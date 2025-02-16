!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CAL_ARPSDERIVED            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cal_arpsderived(imove,umove,vmove,umove_readin,vmove_readin, &
                  presaxis_no,nx,ny,nz,zp,zpc,x,y,mapfct,               &
                  ubar,vbar,wbar,qvbar,ptbar,pbar,                      &
                  uprt,vprt,wprt,qvprt,ptprt,pprt,                      &
                  u,v,w,qv,pt,algpzc,tz,t700,zlevel,psl,td,thet,        &
                  capeplt,cinplt,capsplt,liplt,brnplt,ctcplt,           &
                  li,cape,cin,capst,ct,                                 &
                  heliplt,uhplt,brnuplt,srlfplt,srmfplt,blcoplt,strmplt,&
                  uhmnhgt,uhmxhgt,                                      &
                  ustrm,vstrm,srlfl,srmfl,heli,uh,brn,brnu,blcon,       &
                  ibgnl,iendl,jbgnl,jendl,kbgn,kend,                    &
                  hterain,ztmax,ztmin,u1,v1,                            &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,         &
                  istatus)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: imove, presaxis_no
  REAL,    INTENT(IN)  :: umove,vmove,umove_readin,vmove_readin
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: zp(nx,ny,nz), zpc(nx,ny,nz)
  REAL,    INTENT(IN)  :: x(nx), y(ny)
  REAL,    INTENT(OUT) :: mapfct(nx,ny)

  REAL, DIMENSION(nx,ny,nz), INTENT(IN) :: ubar,vbar,wbar,qvbar,ptbar,pbar
  REAL, DIMENSION(nx,ny,nz), INTENT(IN) :: uprt,vprt,wprt,qvprt,ptprt,pprt
  REAL, DIMENSION(nx,ny,nz), INTENT(OUT) :: u,v,w,qv,pt,algpzc,tz,td
  REAL, DIMENSION(nx,ny),    INTENT(OUT) :: t700,psl,thet
  REAL,                      INTENT(OUT) :: zlevel

  INTEGER, INTENT(IN)  :: capeplt,cinplt,capsplt,liplt,ctcplt
  INTEGER, INTENT(IN)  :: heliplt,uhplt,brnplt,brnuplt
  INTEGER, INTENT(IN)  :: srlfplt,srmfplt,blcoplt,strmplt

  REAL, INTENT(IN) :: uhmnhgt,uhmxhgt
  REAL, DIMENSION(nx,ny), INTENT(OUT) :: li,cape,cin,capst,ct
  REAL, DIMENSION(nx,ny), INTENT(OUT) :: ustrm,vstrm,srlfl,srmfl
  REAL, DIMENSION(nx,ny), INTENT(OUT) :: heli,uh,brn,brnu,blcon

  INTEGER, INTENT(IN)  :: ibgnl,iendl,jbgnl,jendl,kbgn,kend
  REAL,    INTENT(OUT) :: hterain(nx,ny),ztmax,ztmin

  REAL, DIMENSION(nx,ny),    INTENT(INOUT) :: u1,v1
  REAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: tem1,tem2,tem3,tem4,tem5, &
                                              tem6,tem7,tem8,tem9
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: lcl(:,:)    ! The lifting condensation level
  REAL, ALLOCATABLE :: lfc(:,:)    ! Level of free convection
  REAL, ALLOCATABLE :: el(:,:)     ! Equilibrium level
  REAL, ALLOCATABLE :: twdf(:,:)   ! Max wet-bulb potential temperature difference
  REAL, ALLOCATABLE :: mbe(:,:)    ! Moist Convective Potential Energy
                                   !   (MCAPE, includes water loading)

  REAL, ALLOCATABLE :: wrk1(:),wrk2(:),wrk3(:),wrk4(:),wrk5(:),wrk6(:)
  REAL, ALLOCATABLE :: wrk7(:),wrk8(:),wrk9(:),wrk10(:),wrk11(:),wrk12(:)

  REAL, ALLOCATABLE :: shr37(:,:)  ! 7km - 3km wind shear

  REAL, PARAMETER :: gamma=6.5,              & ! 6.5 K/km
                     ex1=0.1903643,          & ! R*gamma/g
                     ex2=5.2558774,          & ! g/R/gamma
                     mbtopa=100.

  REAL    :: uadd, vadd
  REAL    :: p00, t00

  INTEGER :: i,j,k

  REAL    :: f_cputime,cpu0,cpu1,cpu2
  DOUBLE PRECISION :: f_walltime,second1, second2

  INCLUDE 'grid.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
!
!-----------------------------------------------------------------------
!
!  Add the domain translation/grid movement speed back to the wind fields,
!  so that the wind is the Earth relative.
!
!-----------------------------------------------------------------------
!
  IF( imove == 1 ) THEN
    uadd = umove_readin
    vadd = vmove_readin
  ELSEIF( imove == 2 ) THEN
    uadd = umove
    vadd = vmove
  ELSE
    uadd = 0.0
    vadd = 0.0
  END IF

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        u(i,j,k)=uprt(i,j,k)+ubar(i,j,k)+uadd
        v(i,j,k)=vprt(i,j,k)+vbar(i,j,k)+vadd
        w(i,j,k)=wprt(i,j,k)+wbar(i,j,k)
        qv(i,j,k)=MAX(0.0,qvprt(i,j,k)+qvbar(i,j,k))
        pt(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
        tem1(i,j,k)=0.0
        tem2(i,j,k)=0.0
        tem3(i,j,k)=0.0
        tem4(i,j,k)=0.0
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Set negative logrithm pressure (Pascal) coordinates
!  for scalar and w grid points
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        algpzc(i,j,k) = -ALOG(pbar(i,j,k)+pprt(i,j,k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate temperature (K) at ARPS grid points and 700mb
!  pressure level
!
!-----------------------------------------------------------------------
!
  CALL temper (nx,ny,nz,pt, pprt ,pbar,tz)
  zlevel=-ALOG(70000.0)
  CALL hintrp(nx,ny,nz,tz,algpzc,zlevel,t700)
!
!----------------------------------------------------------------------
!
!  Calculate sea level pressure (mb)
!  Reduction method: Benjamin and Miller: 1990, MWR, vol.118, No.10,
!                   Page: 2100-2101
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      p00 = 0.01*(pbar(i,j,2)+pprt(i,j,2))
      IF(p00 <= 700.0) THEN
        t00=tz(i,j,2)
      ELSE
        t00 = t700(i,j)*(p00/700.0)**ex1
      END IF
      psl(i,j) = p00*((t00+gamma*zpc(i,j,2))/t00)**ex2
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate dew-point temperature td using enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
        tem3(i,j,k) = ptbar(i,j,k)+ptprt(i,j,k)
      END DO
    END DO
  END DO

  CALL temper(nx,ny,nz,tem3, pprt ,pbar,tem2)
  CALL getdew(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, tem1, tem2, qv, td)
!
!-----------------------------------------------------------------------
!
!  Calculate CAPE and CIN
!  change qv to mixing ratio
!
!-----------------------------------------------------------------------
!
  CALL edgfill(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1)

!    IF (mp_opt > 0) THEN    ! If we know the read-in data is valid
!                            ! there is no need to call message passing
!      CALL acc_interrupt(mp_acct)
!      CALL mpsendrecv2dew(u,nx,ny,nz,1,1,1,tem4)
!      CALL mpsendrecv2dns(u,nx,ny,nz,1,1,1,tem5)
!
!      CALL mpsendrecv2dew(v,nx,ny,nz,1,1,2,tem4)
!      CALL mpsendrecv2dns(v,nx,ny,nz,1,1,2,tem5)
!
!      CALL mpsendrecv2dew(qv,nx,ny,nz,1,1,0,tem4)
!      CALL mpsendrecv2dns(qv,nx,ny,nz,1,1,0,tem5)
!
!      CALL mpsendrecv2dew(pprt,nx,ny,nz,1,1,0,tem4)
!      CALL mpsendrecv2dns(pprt,nx,ny,nz,1,1,0,tem5)
!
!      CALL mpsendrecv2dew(pbar,nx,ny,nz,1,1,0,tem4)
!      CALL mpsendrecv2dns(pbar,nx,ny,nz,1,1,0,tem5)
!    END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem6(i,j,k)=0.5*(u(i,j,k)+u(i+1,j,k))
        tem7(i,j,k)=0.5*(v(i,j,k)+v(i,j+1,k))
        tem8(i,j,k)=qv(i,j,k)/(1.-qv(i,j,k))
        tem9(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
      END DO
    END DO
  END DO

  CALL edgfill(tem6,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(tem7,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(tem8,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(tem9,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(tz,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
  CALL edgfill(td,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

  ALLOCATE(wrk1 (nz),wrk2 (nz),wrk3 (nz),STAT=istatus)
  wrk1 =0.0; wrk2 = 0.0; wrk3 = 0.0

  IF (capeplt /= 0 .OR. cinplt /= 0 .OR. capsplt /= 0 .OR.            &
        liplt /= 0 .OR. brnplt /= 0) THEN

    second1= f_walltime()
    cpu1 = f_cputime()

    ALLOCATE(lcl  (nx,ny),STAT=istatus)
    ALLOCATE(lfc  (nx,ny),STAT=istatus)
    ALLOCATE(el   (nx,ny),STAT=istatus)
    ALLOCATE(twdf (nx,ny),STAT=istatus)
    ALLOCATE(mbe  (nx,ny),STAT=istatus)

    lcl  =0.0
    lfc  =0.0
    el   =0.0
    twdf =0.0
    mbe  =0.0

    ALLOCATE(wrk4 (nz),wrk5 (nz),wrk6 (nz),STAT=istatus)
    ALLOCATE(wrk7 (nz),wrk8 (nz),wrk9 (nz),STAT=istatus)
    ALLOCATE(wrk10(nz),wrk11(nz),wrk12(nz),STAT=istatus)
    wrk4 =0.0; wrk5 = 0.0; wrk6 = 0.0
    wrk7 =0.0; wrk8 = 0.0; wrk9 = 0.0
    wrk10=0.0; wrk11 = 0.0; wrk12 = 0.0

    CALL arps_be(nx,ny,nz,                                            &
                 tem9,zpc,tz,tem8,                                    &
                 lcl,lfc,el,twdf,li,cape,mbe,cin,capst,               &
                 wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,                       &
                 wrk7,wrk8,wrk9,wrk10,wrk11,wrk12,tem5)

    DEALLOCATE(lcl,lfc,el,twdf,mbe)
    DEALLOCATE(wrk4,wrk5,wrk6)
    DEALLOCATE(wrk7,wrk8,wrk9)
    DEALLOCATE(wrk10,wrk11,wrk12)

    cpu2 = f_cputime()
    second2 = f_walltime()

!    print*,'!!!!  total clock time for arps_be:',second2-second1
!    print*,'!!!!  total cpu time for arps_be  :', cpu2-cpu1
!    print *, ' back from arps_be'
  END IF

!-----------------------------------------------------------------------
!
!  Calculate the convective temperature (celsius)
!
!-----------------------------------------------------------------------

  IF(ctcplt /= 0) THEN
    CALL arps_ct(ct,nx,ny,nz,tem9,tz,td,tem1,wrk1,wrk2,wrk3 )
  END IF

  DEALLOCATE(wrk1,wrk2,wrk3)

!
!-----------------------------------------------------------------------
!
!  Store k=2 theta-e and dewpt
!  level 1 and 2 respectively.
!
!-----------------------------------------------------------------------
!
  CALL pt2pte(nx,ny,1, 1,nx-1,1,ny-1,1,1,                               &
              tem9(1,1,2),tem3(1,1,2),qv(1,1,2),thet)
!
!-----------------------------------------------------------------------
!
!  Calculate helicity and other shear-related paramaters
!
!-----------------------------------------------------------------------
!
  IF (heliplt /= 0 .OR.  brnplt /= 0 .OR. brnuplt /= 0 .OR.           &
      srlfplt /= 0 .OR. srmfplt /= 0 .OR. blcoplt /= 0 .OR.           &
      strmplt /= 0 ) THEN

    ALLOCATE(shr37(nx,ny),STAT=istatus)
    shr37=0.0

    CALL xytomf(nx,ny,x,y,mapfct)
    CALL calcshr(nx,ny,nz,x,y,zp,mapfct,                              &
                 tem9,tz,tem6,tem7,cape,                              &
                 shr37,ustrm,vstrm,srlfl,srmfl,heli,brn,brnu,blcon,   &
                 tem1,u1,v1,tem2,tem5)

    DEALLOCATE(shr37)

  END IF

  IF (uhplt /= 0 ) THEN

    CALL xytomf(nx,ny,x,y,mapfct)
    WRITE(6,'(1x,a,2(F12.2,a))') 'Calculating UH from ',uhmnhgt,' to ',uhmxhgt,' m AGL'
    CALL calc_uh(nx,ny,nz,zp,mapfct,dx,dy,uhmnhgt,uhmxhgt,            &
                 tem6,tem7,w,uh,tem1,tem2)

  END IF
!
!----------------------------------------------------------------------
!
!    Calculate the coordinate value of pres_z to pres_val
!    if presaxis_no >0
!
!----------------------------------------------------------------------
!
  IF( presaxis_no > 0 )                                               &
      CALL interp_p (pbar,zpc,ibgnl,iendl,nx,jbgnl,jendl,ny,kbgn,kend,nz)

!
!-----------------------------------------------------------------------
!
!  Set terrain data
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      hterain(i,j)=zp(i,j,2)       ! Give zp(i,j,2) to hterain(i,j)
    END DO
  END DO

  ztmin=hterain(1,1)
  ztmax=ztmin

  DO j=1,ny
    DO i=1,nx
      ztmax=MAX(ztmax,hterain(i,j))
      ztmin=MIN(ztmin,hterain(i,j))
    END DO
  END DO

  CALL mpmax0(ztmax,ztmin)

  RETURN
END SUBROUTINE cal_arpsderived
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE CALC_UH                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  Calculates updraft helicity (UH) to detect rotating updrafts.
!  Formula follows Kain et al, 2008, Wea. and Forecasting, 931-952,
!  but this version has controls for the limits of integration
!  uhminhgt to uhmxhgt, in m AGL.  Kain et al used 2000 to 5000 m.
!  Units of UH are m^2/s^2.
!
!  Note here that us and vs are at ARPS scalar points.
!  w is at w-point (scalar pt in horiz, staggered vertical)
!
!  Keith Brewster, CAPS/Univ. of Oklahoma
!  March, 2010
!
SUBROUTINE calc_uh(nx,ny,nz,zp,mapfct,dx,dy,uhmnhgt,uhmxhgt,        &
                   us,vs,w,uh,tem1,tem2)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL, INTENT(IN)  :: zp(nx,ny,nz)
  REAL, INTENT(IN)  :: mapfct(nx,ny)
  REAL, INTENT(IN)  :: dx,dy
  REAL, INTENT(IN)  :: uhmnhgt,uhmxhgt
  REAL, INTENT(IN)  :: us(nx,ny,nz)
  REAL, INTENT(IN)  :: vs(nx,ny,nz)
  REAL, INTENT(IN)  :: w(nx,ny,nz)
  REAL, INTENT(OUT) :: uh(nx,ny)
  REAL, INTENT(OUT) :: tem1(nx,ny,nz)
  REAL, INTENT(OUT) :: tem2(nx,ny,nz)
!
! Misc local variables
!
  INTEGER :: i,j,k,kbot,ktop
  REAL    :: twodx,twody,wgtlw,sum,wmean,wsum,wavg
  REAL    :: helbot,heltop,wbot,wtop
  REAL    :: zbot,ztop
!
! Initialize arrays
!
  uh=0.0
  tem1=0.0
!
! Calculate vertical component of helicity at scalar points
!   us: u at scalar points
!   vs: v at scalar points
!
  twodx=2.0*dx
  twody=2.0*dy
  DO k=2,nz-2
    DO j=2,ny-1
      DO i=2,nx-1
        wavg=0.5*(w(i,j,k)+w(i,j,k+1))
        tem1(i,j,k)=wavg *                                      &
            ((vs(i+1,j,k)-vs(i-1,j,k))/(twodx*mapfct(i,j))  -   &
             (us(i,j+1,k)-us(i,j-1,k))/(twody*mapfct(i,j)))
        tem2(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
!
! Integrate over depth uhminhgt to uhmxhgt AGL
!
  WRITE(6,'(1x,a,2(f12.2,a))') &
        'Calculating UH from ',uhmnhgt,' to ',uhmxhgt,' m AGL'
  DO j=2,ny-2
    DO i=2,nx-2
      zbot=zp(i,j,2)+uhmnhgt
      ztop=zp(i,j,2)+uhmxhgt
!
! Find wbar, weighted-mean vertical velocity in column
! Find w at uhmnhgt AGL (bottom)
!
      DO k=2,nz-3
        IF(zp(i,j,k) > zbot) EXIT
      END DO
      kbot=k
      wgtlw=(zp(i,j,kbot)-zbot)/(zp(i,j,kbot)-zp(i,j,kbot-1))
      wbot=(wgtlw*w(i,j,kbot-1))+((1.-wgtlw)*w(i,j,kbot))
!
! Find w at uhmxhgt AGL (top)
!
      DO k=2,nz-3
        IF(zp(i,j,k) > ztop) EXIT
      END DO
      ktop=k
      wgtlw=(zp(i,j,ktop)-ztop)/(zp(i,j,ktop)-zp(i,j,ktop-1))
      wtop=(wgtlw*w(i,j,ktop-1))+((1.-wgtlw)*w(i,j,ktop))
!
! First part, uhmnhgt to kbot
!
      wsum=0.5*(w(i,j,kbot)+wbot)*(zp(i,j,kbot)-zbot)
!
! Integrate up through column
!
      DO k=(kbot+1),(ktop-1)
        wsum=wsum+0.5*(w(i,j,k)+w(i,j,k-1))*(zp(i,j,k)-zp(i,j,k-1))
      END DO
!
! Last part, ktop-1 to uhmxhgt
!
      wsum=wsum+0.5*(wtop+w(i,j,ktop-1))*(ztop-zp(i,j,ktop-1))
      wmean=wsum/(uhmxhgt-uhmnhgt)

      IF(wmean > 0.) THEN    ! column updraft, not downdraft
!
! Find helicity at uhmnhgt AGL (bottom)
!
        DO k=2,nz-3
          IF(tem2(i,j,k) > zbot) EXIT
        END DO
        kbot=k
        wgtlw=(tem2(i,j,kbot)-zbot)/(tem2(i,j,kbot)-tem2(i,j,kbot-1))
        helbot=(wgtlw*tem1(i,j,kbot-1))+((1.-wgtlw)*tem1(i,j,kbot))
!
! Find helicity at uhmxhgt AGL (top)
!
        DO k=2,nz-3
          IF(tem2(i,j,k) > ztop) EXIT
        END DO
        ktop=k
        wgtlw=(tem2(i,j,ktop)-ztop)/(tem2(i,j,ktop)-tem2(i,j,ktop-1))
        heltop=(wgtlw*tem1(i,j,ktop-1))+((1.-wgtlw)*tem1(i,j,ktop))
!
! First part, uhmnhgt to kbot
!
        sum=0.5*(tem1(i,j,kbot)+helbot)*(tem2(i,j,kbot)-zbot)
!
! Integrate up through column
!
        DO k=(kbot+1),(ktop-1)
          sum=sum+0.5*(tem1(i,j,k)+tem1(i,j,k-1))*(tem2(i,j,k)-tem2(i,j,k-1))
        END DO
!
  ! Last part, ktop-1 to uhmxhgt
!
        uh(i,j)=sum+0.5*(heltop+tem1(i,j,ktop-1))*(ztop-tem2(i,j,ktop-1))
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE calc_uh
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SET_COORD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE set_coord(iorig,xgrdorg,ygrdorg,xorig,yorig,                 &
                     nx,ny,nz,nzsoil,nzmax,                             &
                     x,y,z,zp,zpsoil,xc,yc,zc,zpc,zpsoilc,              &
                     dxinv,dyinv,dzinv,dxkm,dykm,dzkm,dzsoilcm,         &
                     xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend,   &
                     iskip,jskip,                                       &
                     ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend,   &
                     xpbgn,xpend,ypbgn,ypend,ibgnl,iendl,jbgnl,jendl,   &
                     istride,jstride,kstride,ist,jst,kst,ibgnlv,jbgnlv, &
                     xr, yr, x1,y1,x2,y2, zr,z1,z2, zsoilr,zsoil1,zsoil2,&
                     zmin,zmax,zsoilmin,zsoilmax,dbglvl,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE: Set plotting coordinates and domain index boundes
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iorig
  REAL,    INTENT(IN) :: xgrdorg, ygrdorg
  REAL,    INTENT(IN) :: xorig, yorig
  INTEGER, INTENT(IN) :: nx,ny,nz,nzsoil,nzmax
  REAL,    INTENT(IN) :: x(nx),y(ny),z(nz)
  REAL,    INTENT(IN) :: zp(nx,ny,nz), zpsoil(nx,ny,nzsoil)
  REAL,    INTENT(OUT) :: xc(nx,ny,nzmax),yc(nx,ny,nzmax)
  REAL,    INTENT(OUT) :: zc(nx,ny,nz),zpc(nx,ny,nz),zpsoilc(nx,ny,nzsoil)

  REAL,    INTENT(OUT) :: dxinv, dyinv, dzinv
  REAL,    INTENT(OUT) :: dxkm,dykm, dzkm,dzsoilcm

  INTEGER, INTENT(IN)    :: iskip, jskip
  REAL,    INTENT(INOUT) :: xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend
  INTEGER, INTENT(OUT)   :: ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend
  INTEGER, INTENT(OUT)   :: xpbgn,xpend,ypbgn,ypend,ibgnl,iendl,jbgnl,jendl

  INTEGER, INTENT(IN)    :: istride,jstride, kstride   ! For vector plot only
  INTEGER, INTENT(OUT)   :: ist, jst, kst, ibgnlv, jbgnlv

  REAL,    INTENT(OUT)   :: xr, yr, x1,y1,x2,y2, zr,z1,z2, zsoilr,zsoil1,zsoil2
  REAL,    INTENT(OUT)   :: zmin,zmax,zsoilmin,zsoilmax

  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  REAL    :: dx, dy, dz
  REAL    :: xorold, yorold
  INTEGER :: procbgn

  INTEGER :: i, j, k
  INTEGER :: idiff, jdiff

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
!
!-----------------------------------------------------------------------
!
!  Set coordinates at the grid volume center
!
!-----------------------------------------------------------------------
!
  DO k=1,nzmax
    DO j=1,ny
      DO i=1,nx-1
        xc(i,j,k) = (x(i)+x(i+1))*0.5 *0.001
      END DO
    END DO
  END DO

  DO k=1,nzmax
    DO j=1,ny-1
      DO i=1,nx
        yc(i,j,k) = (y(j)+y(j+1))*0.5 *0.001
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx
        zc(i,j,k) = (z(k)+z(k+1))*0.5 *0.001
        zpc(i,j,k)= (zp(i,j,k)+zp(i,j,k+1))*0.5 *0.001  ! km
      END DO
    END DO
  END DO

!
! 05/30/2002 Zuwen He
!
! zpsoilc is the vertical distance of the soil model
! from the terrain surface in meters.
!

  DO k=1,nzsoil
    DO j=1,ny
      DO i=1,nx
        zpsoilc(i,j,k)=(zpsoil(i,j,k)-zpsoil(i,j,1))*100.
                      ! cm for soil model, Zuwen
      END DO
    END DO
  END DO
!
! END Zuwen
!
!
!-----------------------------------------------------------------------
!
!  Modify the coordinate arrays to make certain that the origin is
!  the same as specified above.
!
!-----------------------------------------------------------------------
!
    dx = x(2)-x(1) ! in meters
    dy = y(2)-y(1) ! in meters
    dz = z(2)-z(1) ! in meters
    dxinv = 1.0/(x(2)-x(1))
    dyinv = 1.0/(y(2)-y(1))
    dzinv = 1.0/(z(2)-z(1))

    IF( iorig == 1) THEN

      xc = xc + xgrdorg * 0.001
      yc = yc + ygrdorg * 0.001

    ELSE IF( iorig == 2) THEN

      xorold = (xc(1,2,2)+xc(2,2,2))/2 + xgrdorg * 0.001
      yorold = (yc(2,1,2)+yc(2,2,2))/2 + ygrdorg * 0.001

      dxkm = xc(3,2,2)-xc(2,2,2) ! dx in km
      dykm = yc(2,3,2)-yc(2,2,2) ! dy in km

      WRITE(6,'(/1x,4(a,f8.3),a/)')                                     &
          'Model grid origin (',xorold,',',yorold,') was reset to ('    &
                          ,xorig,',',yorig,').'

      IF( ABS(xorig-xorold) > 0.01*dxkm .OR.                            &
          ABS(yorig-yorold) > 0.01*dykm ) THEN

        xc = xc - xorold + xorig
        yc = yc - yorold + yorig

      END IF

    ELSE IF( iorig == 3) THEN

      xc = xc + xorig
      yc = yc + yorig

    END IF  ! iorig

!-----------------------------------------------------------------------
!
!  Get index bounds of the domain to be plotted
!
!-----------------------------------------------------------------------
!
    CALL indxbnds(xc,yc,zpc,zpsoilc,nx,ny,nz,nzsoil,                    &
                  xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend,      &
                  ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend)

!ibgn, iend, jbgn, jend, kbgn, kend etc. computed above are over global domain
!Introduced ibgnl,iendl, jbgnl,jendl below for message passing
!version, which are relative to the local domain.   - WYH

    xpbgn = (ibgn-2)/(nx-3) + 1
    xpend = (iend-2)/(nx-3) + 1
    ypbgn = (jbgn-2)/(ny-3) + 1
    ypend = (jend-2)/(ny-3) + 1

    IF (loc_x == xpbgn) THEN
      ibgnl = MOD((ibgn-2),(nx-3)) + 2   ! convert to local index
    ELSE IF (loc_x < xpbgn .OR. loc_x > xpend) THEN
      ibgnl = 2          ! can be any value
    ELSE
      ibgnl = 2          ! should be 2, 1 for overlay
    END IF

    IF (loc_x == xpend) THEN
      iendl = MOD((iend-2),(nx-3)) + 2
    ELSE IF (loc_x > xpend .OR. loc_x < xpbgn) THEN
      iendl = nx-2       ! can be any value
    ELSE
      iendl = nx - 1     ! should be nx-2, nx-1 for overlay
    END IF

    IF (loc_y == ypbgn) THEN
      jbgnl = MOD((jbgn-2),(ny-3)) + 2
    ELSE IF (loc_y < ypbgn .OR. loc_y > ypend) THEN
      jbgnl = 2
    ELSE
      jbgnl = 2
    END IF

    IF (loc_y == ypend) THEN
      jendl = MOD((jend-2),(ny-3)) + 2
    ELSE IF (loc_y > ypend .OR. loc_y < ypbgn) THEN
      jendl = ny-2
    ELSE
      jendl = ny - 1
    END IF

    IF (dbglvl > 1) THEN
      IF (myproc == 0) WRITE(6,'(1x,a,4(I4,a))')                        &
      '--- The global range to be plotted is ',ibgn,' -',iend,', ',jbgn,' -',jend,'.'
      WRITE(6,'(1x,a,5(I4,a))') '--- Range of processor -',myproc,      &
        '       is ',ibgnl,' -',iendl,', ',jbgnl, ' -',jendl,'.'
    END IF
! END  --  WYH

    IF( istride /= 0 ) THEN
      ist = istride
    ELSE IF ( (iend-ibgn+1) > 110) THEN
      ist = 4
    ELSE IF( (iend-ibgn+1) > 30) THEN
      ist = 2
    ELSE
      ist = 1
    END IF

    IF( jstride /= 0 ) THEN
      jst = jstride
    ELSE IF ( (jend-jbgn+1) > 110) THEN
      jst = 4
    ELSE IF( (jend-jbgn+1) > 30) THEN
      jst = 2
    ELSE
      jst = 1
    END IF

    ibgnlv = (loc_x-1)*(nx-3)+ibgnl  ! hold global index for ibgnl temporarily
    jbgnlv = (loc_y-1)*(ny-3)+jbgnl  ! hold global index for jbgnl
    IF(MOD(ibgnlv-ibgn,ist) /= 0)   THEN
       ibgnlv = ibgnl + (ist- MOD(ibgnlv - ibgn,ist)) ! for vector plot
    ELSE
       ibgnlv = ibgnl
    END IF
    IF(MOD(jbgnlv-jbgn,jst) /= 0)   THEN
       jbgnlv = jbgnl + (jst- MOD(jbgnlv - jbgn,jst)) ! for vector plot
    ELSE
       jbgnlv = jbgnl
    END IF

    IF( kstride /= 0 ) THEN
      kst = kstride
    ELSE IF ( (kend-kbgn+1) > 90 ) THEN
      kst = 3
    ELSE IF( (kend-kbgn+1) > 45) THEN
      kst = 2
    ELSE
      kst = 1
    END IF

!-----------------------------------------------------------------------
!
! Adjustment of ibgnl and jbgnl based on iskip & jskip.
!
! NOTE: Vector plot is not affected because it uses ibgnlv, jbgnlv above.
!
!-----------------------------------------------------------------------

!  IF (iskip > 0 .AND. loc_x > xpbgn .AND. loc_x <= xpend) THEN
!    idiff = MOD( ((nx-3)*(loc_x-1)+ibgnl-ibgn), (iskip+1) )
!    IF (idiff /= 0) THEN
!      ibgnl = ibgnl + (iskip + 1 - idiff)
!      IF (dbglvl > 1) WRITE(6,'(1x,2(a,I4),a,I2)')                    &
!      '%%% Processor - ',myproc,' reset ibgnl = ',ibgnl,' based on iskip = ',iskip
!    END IF
!  END IF
!
!  IF (jskip > 0 .AND. loc_y > ypbgn .AND. loc_y <= ypend) THEN
!    jdiff = MOD( ((ny-3)*(loc_y-1)+jbgnl-jbgn), (jskip+1) )
!    IF (jdiff /= 0) THEN
!      jbgnl = jbgnl + (jskip + 1 - jdiff)
!      IF (dbglvl > 1) WRITE(6,'(1x,2(a,I4),a,I2)')                    &
!      '%%% Processor - ',myproc,' reset jbgnl = ',jbgnl,' based on jskip = ',jskip
!    END IF
!  END IF

!--------------- End of adjustment -------------------------------------

  dxkm = xc(3,2,2)-xc(2,2,2)
  dykm = yc(2,3,2)-yc(2,2,2)
  dzkm = zc(2,2,3)-zc(2,2,2)
  dzsoilcm = ABS(zpsoilc(1,1,2)-zpsoilc(1,1,1))/2.

  xr = (iend-ibgn)*dxkm + dxkm
  yr = (jend-jbgn)*dykm + dykm

  procbgn = (ypbgn-1)*nproc_x+xpbgn-1  ! processor rank which holds (ibgn,jbgn)
  IF (myproc == procbgn) THEN
    x1 = xc(ibgnl,2,2)-dxkm/2
    y1 = yc(2,jbgnl,2)-dykm/2
  END IF
  CALL mpbcastr(x1,procbgn)
  CALL mpbcastr(y1,procbgn)

  x2 = x1 + xr
  y2 = y1 + yr

  IF( zbgn /= zend ) THEN
    zmin=zbgn     ! these values are already global
    zmax=zend
  ELSE

    IF (loc_x >= xpbgn .AND. loc_x <= xpend .AND.                     &
        loc_y >= ypbgn .AND. loc_y <= ypend) THEN

      zmax=zp(ibgnl,jbgnl,kend+1)
      zmin=zp(ibgnl,jbgnl,kbgn)

      DO j=jbgnl,jendl
        DO i=ibgnl,iendl
          zmax=AMAX1(zmax,zp(i,j,kend+1))
          zmin=AMIN1(zmin,zp(i,j,kbgn)  )
        END DO
      END DO

      IF( kbgn == 2 ) zmin=z(kbgn)

      zmin=zmin / 1000.
      zmax=zmax / 1000.

    ELSE         ! we do not care their values
      zmin =  1.0E6    ! should be large enough
      zmax = -1.0E6    ! should be small enough
    END IF

    CALL mpmax0(zmax,zmin)    ! get global values

  END IF

  zr = zmax-zmin
  z1 = zmin
  z2 = zmax

!
! 05/31/2002 Zuwen He
!
! For soil model, the zsoilmax and zsoilmin are the
! maximum and minimum of zpsoilc.
!

  IF( zsoilbgn /= zsoilend ) THEN
    zsoilmax=zsoilbgn
    zsoilmin=zsoilend
  ELSE
    IF (loc_x >= xpbgn .AND. loc_x <= xpend .AND.                     &
        loc_y >= ypbgn .AND. loc_y <= ypend) THEN

      zsoilmax=zpsoilc(ibgnl,jbgnl,ksoilbgn)
      zsoilmin=zpsoilc(ibgnl,jbgnl,ksoilend)

      DO j=jbgnl,jendl
        DO i=ibgnl,iendl
          zsoilmax=AMAX1(zsoilmax,zpsoilc(i,j,ksoilbgn))
          zsoilmin=AMIN1(zsoilmin,zpsoilc(i,j,ksoilend))
        END DO
      END DO
    ELSE         ! we do not care their values
      zsoilmin =  1.0E6    ! should be large enough
      zsoilmax = -1.0E6    ! should be small enough
    END IF

    CALL mpmax0(zsoilmax,zsoilmin)

  END IF

  zsoilr = max(zsoilmax-zsoilmin,0.0001)
  zsoil1 = zsoilmin
  zsoil2 = zsoilmax

  RETURN
END SUBROUTINE set_coord

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READ_STATION               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE read_station(infile,mxstalo,latsta,lonsta,                   &
           nstatyp,nstapro,nsta,sname,state,sitena,nelev)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will read external staion information.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (6/1/97)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=*) :: infile
  INTEGER :: nsta,nstapro(mxstalo),nstatyp(mxstalo)
  REAL :: latsta(mxstalo), lonsta(mxstalo)
  CHARACTER (LEN=5) :: sname(mxstalo)
  CHARACTER (LEN=2) :: state(mxstalo)
  CHARACTER (LEN=20) :: sitena(mxstalo)
  INTEGER :: nelev(mxstalo)
  CHARACTER (LEN=132) :: line

  OPEN(1,IOSTAT=ios,FILE=infile,STATUS='old',                           &
        FORM='formatted')
  IF(ios /= 0) THEN     ! error during read
    istatus = -1
    WRITE(6,650) infile
    650       FORMAT(' +++ ERROR opening: ',a70,' +++')
    WRITE(6,651) ios
    651       FORMAT('     IOS code = ',i5)
    RETURN
  END IF

  nsta = 0

! Read only lines that begin with A-Z, a-z, or 0-9 -- treat the rest as
! comments

  DO i=1,mxstalo
    READ(1,'(A)',END=999) line
    IF ( ( line(:1) >= 'A' .AND. line(:1) <= 'Z' ) .OR.                 &
           ( line(:1) >= 'a' .AND. line(:1) <= 'z' ) .OR.               &
           ( line(:1) >= '0' .AND. line(:1) <= '9' ) ) THEN
      nsta = nsta + 1
      READ(line,101,ERR=999)sname(nsta),state(nsta),sitena(nsta),       &
          latsta(nsta),lonsta(nsta),nelev(nsta),nstatyp(nsta),          &
          nstapro(nsta)
    END IF
    101   FORMAT(a5,2X,a2,1X,a20,4X,f9.4,f9.4,i5,2X,i2,i1)
    102   FORMAT(a5,2X,a2,1X,a20,4X,f9.4,f9.4,i5,2X,i2,1X,i1)
  END DO
  999   CONTINUE
  CLOSE(1)

  RETURN
END SUBROUTINE read_station
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTERP_P                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE interp_p (pbar,zpc,ibgn,iend,nx,jbgn,jend,ny,                &
                     kbgn,kend,nz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will interpolate pressure for draw pressure bar.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (6/1/97)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!

  REAL :: pbar(nx,ny,nz), zpc(nx,ny,nz)
  INTEGER :: presaxis_no
  REAL :: pres_val(20), pres_z(20)
  COMMON /pressbar_par/presaxis_no,pres_val,pres_z
  REAL :: tmp
  REAL :: pres_val1(20)
  REAL :: pz(100), pb(100)

  REAL :: pbarmax, zpcmax
  INTEGER :: ilocs, iloce, jlocs, jloce

!----------------------------------------------------------------------
!
! Include files
!
!----------------------------------------------------------------------
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! estimate a value for pz(kbgn-1) and pb(kbgn-1)
  kbgn1     = kbgn-1
  pb(kbgn1) = 1010.
  pz(kbgn1) = 0.5

! get max pressure

  ilocs = ibgn
  iloce = iend
  jlocs = jbgn
  jloce = jend

  IF (loc_x > 1)       ilocs = ibgn+1  ! ensure no overlap
  IF (loc_x < nproc_x) iloce = iend-1
  IF (loc_y > 1)       jlocs = jbgn+1  ! ensure no overlap
  IF (loc_y < nproc_y) jloce = jend-1

  DO k=kbgn,kend
      pb(k)=0.
      pz(k)=0.
      ini = 0
      inj = 0

      DO j = jlocs, jloce
        DO i = ilocs, iloce
          IF(0.01*pbar(i,j,k) > pb(k)) THEN
            ini = i
            inj = j
            pb(k) = 0.01*pbar(i,j,k)  ! Add by WYH, it was missed originally
          END IF
        END DO
      END DO

      IF(ini /= 0 .AND. inj /= 0) THEN   ! local ini,inj
         pbarmax = 0.01*pbar(ini,inj,k)
      ELSE                     !! missing values
         pbarmax = -1.0
      END IF

      CALL globalpbar(pbarmax,ini,inj,k,zpc,nx,ny,nz,zpcmax)

    IF(ini /= 0 .AND. inj /= 0) THEN   ! global ini, inj
!        pb(k) = 0.01*pbar(ini,inj,k)
!        pz(k) = zpc(ini, inj,k)
         pb(k) = pbarmax
         pz(k) = zpcmax

    ELSE                     !! missing values
        WRITE(6,'(a,i2,a)')                                               &
            'Warning: Missing pressure value on level ',k,                &
            ' Using previous level instead of.'
        pb(k) = pb(k-1)
        pz(k) = pz(k-1)
    END IF
  END DO

  k=0
  DO j=1,presaxis_no
    IF (pres_val(j) > pb(kbgn) ) THEN
      k  = 1
!       print*,'pres_val(j)',pres_val(j)
!       print*, '>= pb(kbgn)',pb(kbgn)
      pres_val1(k) = pb(kbgn)
      pres_z(k) = pz(kbgn)
    END IF
  END DO

  DO j=1,presaxis_no
    DO i=kbgn,kend-1
      tmp = pres_val(j)
      IF(tmp <= pb(i).AND.tmp > pb(i+1))THEN !find pressure interpolate
        k = k+1
        a1 = ALOG(pb(i))-ALOG(tmp)
        a2 = ALOG(pb(i))-ALOG(pb(i+1))
        a3 = pz(i+1)-pz(i)
        pres_z(k) = pz(i)+a3*a1/a2
        pres_val1(k) = pres_val(j)
      END IF
    END DO
  END DO
  k=k+1
  pres_val1(k) = pb(kend)
  pres_z(k) = pz(kend)

  presaxis_no = k
  DO i=1,presaxis_no
    pres_val(i) = pres_val1(i)
  END DO

  RETURN
END SUBROUTINE interp_p
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_MULOVRLAY              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_mulovrlay(var,LEN,num,ovrname,sovrlay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will find out the overlay multiple plots.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (6/1/97)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  INTEGER :: num, sovrlay
  CHARACTER (LEN=*) :: var
  CHARACTER (LEN=*) :: ovrname(10)

  INTEGER :: i,LEN

  sovrlay = 0

  DO i = 1,num
    IF(var(1:LEN) == ovrname(i)(1:LEN) ) sovrlay=1
  END DO

  RETURN
END SUBROUTINE get_mulovrlay

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARPS_CT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE arps_ct (cct,nx,ny,nz,p,t,td,ppccl,wrk1,wrk2,wrk3)

!-----------------------------------------------------------------------
!  Purpose:
!  Calculate the convective temperature (celsius)
!
!  AUTHOR:  Min Zou
!    07/10/1997
!
!  MODIFICATIONS:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    p        pressure (Pascals)
!    t        temperature(degrees Kelvin)
!    td       dew-point temperature (degrees Kelvin)
!
!  OUTPUT:
!
!    cct      convective temperature (Celsius)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  REAL :: p(nx,ny,nz)    ! pressure (Pa)
  REAL :: t(nx,ny,nz)    ! temperature(K)
  REAL :: td(nx,ny,nz)   ! dew-point temperature (K)

  REAL :: cct(nx,ny)     ! convective temperature (C)

!-----------------------------------------------------------------------
!
!  Misc. temporary variables
!
!-----------------------------------------------------------------------

  REAL :: ppccl(nx,ny)   ! pressure (millibars) at the convective
                         !condensation level

  REAL :: wrk1(nz), wrk2(nz), wrk3(nz)

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: pm             ! pressure(millibars) at upper boundary of the
                         ! layer for computing the mean mixing ratio.
  REAL :: mrbar          ! mean mixing ratio (g/kg) in the layer bounded
                         ! by pressures at the p bottom and the pm at the top
!
!-----------------------------------------------------------------------
!
!  Function f_pccl and f_ct and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_pccl, f_ct

!fpp$ expand (f_pccl)
!fpp$ expand (f_ct)
!!dir$ inline always f_pccl,f_ct
!*$*  inline routine (f_pccl,f_ct)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j = 1,ny
    DO i = 1,nx
      pm = p(i,j,2) - 5000.

      ppccl(i,j) = f_pccl(pm,p(i,j,1),t(i,j,1),td(i,j,1),mrbar,nz)

      cct(i,j) = f_ct(mrbar, ppccl(i,j), p(i,j,2))

    END DO
  END DO

  RETURN
END SUBROUTINE arps_ct


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TEMPER                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE temper ( nx,ny,nz,theta, ppert, pbar, t )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Using a version of Poisson's formula, calculate temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Joe Bradley
!    12/05/91
!
!  MODIFICATIONS:
!    Modified by Ming Xue so that arrays are only defined at
!             one time level.
!    6/09/92  Added full documentation and phycst include file for
!             rddcp=Rd/Cp  (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    theta    Potential temperature (degrees Kelvin)
!    ppert    Perturbation pressure (Pascals)
!    pbar     Base state pressure (Pascals)
!
!  OUTPUT:
!
!    t        Temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
!
  REAL :: theta(nx,ny,nz)      ! potential temperature (degrees Kelvin)
  REAL :: ppert(nx,ny,nz)      ! perturbation pressure (Pascals)
  REAL :: pbar (nx,ny,nz)      ! base state pressure (Pascals)
!
  REAL :: t    (nx,ny,nz)      ! temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
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
!  Calculate the temperature using Poisson's formula.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        t(i,j,k) = theta(i,j,k) *                                       &
             (((ppert(i,j,k) + pbar(i,j,k)) / p0) ** rddcp)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE temper

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTEPO                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE intepo(ns,xpos,ypos,vpos,m,n,x,y,var)
!
!  Does bilinear interpolation of variables (nvar variables)
!  to position (xpos,ypos) given a grid of variables "var", their
!  x and y positions of the grid (x,y) and the index of the
!  grid point (i,j) which is to the "lower-left" of xpos and ypos.
!  That means (xpos,ypos) is between (i,j) and (i+1,j+1).
!
!  A saftey feature in case xpos and ypos are outside the grid
!  is built-in.
!
  IMPLICIT NONE
!
!  Arguments, input
!
  INTEGER :: m,n
  REAL :: var(m,n)
  REAL :: x(m,n),y(m,n)
  INTEGER :: ns
  REAL :: xpos(ns),ypos(ns), vpos(ns)
  INTEGER :: i,j
  REAL :: a1,a2
!
!  Arguments output
!
!  Misc Internal Variables
!
  REAL :: dxpos,dypos,dx,dy
  INTEGER :: k
!

  DO k=1,ns
    DO j=1,n-1
      DO i=1,m-1
        IF( (xpos(k) >= x(i,j) .AND. xpos(k) < x(i+1,j+1)) .AND.        &
              (ypos(k) >= y(i,j) .AND. ypos(k) < y(i+1,j+1)) )  THEN

          IF(var(i,j) /= -9999. .AND. var(i+1,j) /= -9999. .AND.        &
                var(i,j+1) /= -9999. .AND. var(i+1,j+1) /= -9999.) THEN
            dx=x(i+1,j)-x(i,j)
            dy=y(i,j+1)-y(i,j)
            dxpos=(xpos(k)-x(i,j))/dx
            dypos=(ypos(k)-y(i,j))/dy

            a1=(1.-dxpos)*var(i,j) + dxpos*var(i+1,j)
            a2=(1.-dxpos)*var(i,j+1) + dxpos*var(i+1,j+1)

            vpos(k) = (1.-dypos)*a1 + dypos*a2
          ELSE
            vpos(k) = -9999.
          END IF
!
        END IF
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE intepo

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UNIGRID                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE unigrid(nx,nz,f,z,fdata,zdata,fprof,zprof)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  E
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       first dimension of f
!    nz       second dimension of f
!    f        2-dimension array of variable
!    z        z coordinate of grid points in physcal space (m)
!    fdata    1-D array defined at levels in zdat to be
!               interpolated to levels defined by zpro.
!    zdata    The height of the input data given by fdat.
!    zprof    The grid level height to which data are interpolated
!
!  OUTPUT:
!
!    fprof    The number of interpolated data levels
!    f        2-dimension array of variable
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,nz
  INTEGER :: i,k
  REAL :: f(nx,nz)
  REAL :: z(nx,nz)
  REAL :: fdata(nz),zdata(nz),fprof(nz),zprof(nz)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=1,nx

    DO k=1,nz
      fdata(k)=f(i,k)
      zdata(k)=z(i,k)
    END DO

    CALL inte1d(fdata,zdata,nz,fprof,zprof,nz)

    DO k=1,nz
      f(i,k)=fprof(k)
      IF(zprof(k) < zdata(1) .OR. zprof(k) > zdata(nz) ) f(i,k)=-9999.0
    END DO

  END DO

  RETURN
END SUBROUTINE unigrid
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FILZERO                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE filzero( a, n)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Fill vector a with zeros.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: St. Paul, Dead Sea Scrolls
!
!  MODIFICATIONS:
!    6/09/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: n
  REAL :: a(n)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=1,n
    a(i)=0.0
  END DO

  RETURN
END SUBROUTINE filzero
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IFILZERO                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ifilzero( ia, n )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Fill vector a with zeros.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: St. Paul, Dead Sea Scrolls
!
!  MODIFICATIONS:
!    6/09/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: n
  INTEGER :: ia(n)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=1,n
    ia(i)=0
  END DO

  RETURN
END SUBROUTINE ifilzero

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SLENGTH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE slength ( string, length )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the length of the non-blank part of a string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  6/09/92 (K. Brewster)
!  Added full documentation and streamlined logic
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   character string to be sized
!
!  INPUT/OUTPUT:
!
!    length   on input, full size of character string
!             on output, true length of string as measured by the
!                        location of the last non-blank character
!
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
  INTEGER :: length
  CHARACTER (LEN=*) :: string
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i = length,1,-1
    IF(string(i:i) /= ' ') EXIT
  END DO
  101   CONTINUE

  length = i

  RETURN
END SUBROUTINE slength
!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_t                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE cal_t ( tem,tz,nx,ny,nz,tob,label,length,units )
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate T value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp, obs_marksz

  INTEGER :: nobs
  COMMON /sfc_obs1/ nobs
  REAL :: latob(mxsfcob),lonob(mxsfcob)
  REAL :: obs1(mxsfcob),obs2(mxsfcob)
  COMMON /sfc_obs2/ latob,lonob,obs1,obs2

  INTEGER :: nx,ny,nz, length
!   real tem(*), tz(*), tob(*)
  REAL :: tem(nx,ny,nz), tz(nx,ny,nz), tob(*)
  CHARACTER (LEN=*) :: units
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j,k, iob
  INTEGER :: ibgn,iend, jbgn, jend, kbgn,kend, isize,jsize,ksize
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ibgn = 1
  iend = nx-1
  jbgn = 1
  jend = ny-1
  kbgn = 1
  kend = nz-1
  isize = (iend-ibgn)+1
  jsize = (jend-jbgn)+1
  ksize = (kend-kbgn)+1

  IF (units(1:1) == 'F' .OR. units(1:1) == 'f') THEN
    DO k= kbgn, kend
      DO j= jbgn, jend
        DO i= ibgn, iend
!       ij = i-ibgn+1 + (j-jbgn)*isize + (k-kbgn)*jsize
!       tem(ii) = 32.0 + 1.8*(tz(ii) - 273.15)
          tem(i,j,k) = 32.0 + 1.8*(tz(i,j,k) - 273.15)
        END DO
      END DO
    END DO
    label = 'T (F)'
    length = 5

    IF(ovrobs == 1 .AND. nobs > 0) THEN
      DO iob=1,nobs
        obs1(iob)=tob(iob)
      END DO
      obsset=1
    END IF

  ELSE                 ! default units is C or c
    DO k= kbgn, kend
      DO j= jbgn, jend
        DO i= ibgn, iend
!       ij = i-ibgn+1 + (j-jbgn)*isize + (k-kbgn)*jsize
!       tem(ii) = tz(ii) - 273.15
          tem(i,j,k) = tz(i,j,k) - 273.15
        END DO
      END DO
    END DO
    label = 'T (C)'
    length = 5

    IF(ovrobs == 1 .AND. nobs > 0) THEN
      DO iob=1,nobs
        obs1(iob)=(tob(iob)-32.)*5./9.
      END DO
      obsset=1
    END IF
  END IF

  RETURN
END SUBROUTINE cal_t

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vh                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate vh value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!  10/17/08 Re-wrote subroutine.  Previous calculation of vh was
!           incorrect.  U and V are not at the same location, so
!           they need to be averaged to the same points first,
!           before calculating the wind speed. (Dan Dawson)
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vh(tem9,u,v,nx,ny,nz,vhunits,label,length,tem8)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz),u(nx,ny,nz), v(nx,ny,nz)
  REAL :: tem8(nx,ny,nz)
  INTEGER :: vhunits, length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j,k, onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ! First average u and v to the scalar points
  onvf = 0
  CALL avgx(u,onvf,nx,ny,nz,2,nx-1,1,ny-1,1,nz-1,tem8)
  CALL avgy(v,onvf,nx,ny,nz,1,nx-1,2,ny-1,1,nz-1,tem9)

  ! Now calculate wind speed
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k) = SQRT(tem8(i,j,k)**2+tem9(i,j,k)**2)
      END DO
    END DO
  END DO

  IF(vhunits == 1) THEN
    label = 'Horiz. wind (m/s)'
    length =17
  ELSE IF( vhunits == 2) THEN
    label = 'Horiz. wind (kts)'
    length = 17
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem9(i,j,k) = tem9(i,j,k)*1.943844
        END DO
      END DO
    END DO
  ELSE IF (vhunits == 3) THEN
    label = 'Horiz. wind (MPH)'
    length = 17
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem9(i,j,k) = tem9(i,j,k)*2.236936
        END DO
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE cal_vh


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_qw                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate qw value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_qw(tem9,qscalar, nx,ny,nz)

  IMPLICIT NONE
  INCLUDE 'globcst.inc'
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k) = 0.0
        IF (P_QC > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QC)
        IF (P_QR > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QR)
        IF (P_QI > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QI)
        IF (P_QS > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QS)
        IF (P_QG > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QG)
        IF (P_QH > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QH)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_qw

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_rh                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate rh value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_rh(tem9,pt, pprt ,pbar,qv,tem1,tem2,nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz), tem1(nx,ny,nz), tem2(nx,ny,nz)
  REAL :: pt(nx,ny,nz), pprt(nx,ny,nz), pbar(nx,ny,nz), qv(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL temper (nx,ny,nz,pt, pprt ,pbar,tem1)
  CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, pbar,tem1,tem2)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem9(i,j,k) = MIN( MAX( qv(i,j,k)/tem2(i,j,k), 0.0), 1.0)
        ! add MIN, MAX as Kevin Thomas' recommendation
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_rh
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_td                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate td value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_td(tem9,td,nx,ny,nz,tdunits,label, length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz), td(nx,ny,nz)
  INTEGER :: length
  CHARACTER (LEN=*) :: label
  CHARACTER (LEN=*) :: tdunits

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (tdunits(1:1) == 'F' .OR. tdunits(1:1) == 'f') THEN
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem9(i,j,k) = 32.0 + 1.8*(td(i,j,k) - 273.15)
        END DO
      END DO
    END DO
    label = 'Td (F)'
    length = 6
  ELSE
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem9(i,j,k) = td(i,j,k) - 273.15
        END DO
      END DO
    END DO
    label = 'Td (C)'
    length = 6
  END IF

  RETURN
END SUBROUTINE cal_td

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_tdobs                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate td observation value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_tdobs(tdob, tdunits)

  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp, obs_marksz

  INTEGER :: nobs
  COMMON /sfc_obs1/ nobs
  REAL :: latob(mxsfcob),lonob(mxsfcob)
  REAL :: obs1(mxsfcob),obs2(mxsfcob)
  COMMON /sfc_obs2/ latob,lonob,obs1,obs2

  REAL :: tdob(mxsfcob)

  INTEGER :: iob
  CHARACTER (LEN=*) :: tdunits
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (tdunits(1:1) == 'F' .OR. tdunits(1:1) == 'f') THEN
    IF(ovrobs == 1 .AND. nobs > 0) THEN
      DO iob=1,nobs
        obs1(iob)=tdob(iob)
      END DO
      obsset=1
    END IF
  ELSE
    IF(ovrobs == 1 .AND. nobs > 0) THEN
      DO iob=1,nobs
        obs1(iob)=(tdob(iob)-32.)*5./9.
      END DO
      obsset=1
    END IF
  END IF

  RETURN
END SUBROUTINE cal_tdobs

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_rfc                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate rfc value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!  Ming Xue (10/16/2001)
!  Now passing in precalculated reflectivity field instead of calculating
!  it inside.
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_rfc(nx, ny, nz, ref, refc)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL, INTENT(IN ) :: ref (nx,ny,nz) ! Reflectivity
  REAL, INTENT(OUT) :: refc(nx,ny,nz) ! Composite reflectivity

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j=1,ny
    DO i=1,nx
      refc(i,j,1)= ref(i,j,1)
      DO k=2,nz-1
        refc(i,j,1) = MAX(refc(i,j,1),ref(i,j,k))
      END DO
    END DO
  END DO

  DO j=1,ny
    DO i=1,nx
      DO k=2,nz-1
        refc(i,j,k) = refc(i,j,1)
      END DO
    END DO
  END DO


  RETURN
END SUBROUTINE cal_rfc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vorp                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vort*10^5 (1/s) value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vorp(tem9,u,v,x,y,nx,ny,nz,tem1)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)
  REAL :: tem1(nx,ny,nz)  ! work array

  INTEGER :: i,j,k

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=2,nz-2
    DO  j=2,ny-2
      DO  i=2,nx-2
        tem9(i,j,k)= 1.0E5*(                                            &
               (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
               (4*(x(i+1)-x(i)))-                                       &
               (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
               (4*(y(j+1)-y(j))) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      tem9(i,j,   1)=tem9(i,j,   2)
      tem9(i,j,nz-1)=tem9(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      tem9(   1,j,k)=tem9(   2,j,k)
      tem9(nx-1,j,k)=tem9(nx-2,j,k)
    END DO
  END DO

  DO  k=1,nz-1
    DO i=1,nx-1
      tem9(i,   1,k)=tem9(i,   2,k)
      tem9(i,ny-1,k)=tem9(i,ny-2,k)
    END DO
  END DO

  IF(mp_opt > 0) THEN
    CALL mpsendrecv2dew(tem9,nx,ny,nz,1,1,0,tem1)
    CALL mpsendrecv2dns(tem9,nx,ny,nz,1,1,0,tem1)
  END IF

  RETURN
END SUBROUTINE cal_vorp

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE cal_div                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate 1000.*Divergence (1/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_div(tem9,u,v,x,y,nx,ny,nz,tem1)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)
  REAL :: tem1(nx,ny,nz)  ! work array

  INTEGER :: i,j,k

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-1
        tem9(i,j,k)=1000.*( (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))         &
                          + (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j)) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      tem9(i,j,   1)=tem9(i,j,   2)
      tem9(i,j,nz-1)=tem9(i,j,nz-2)
    END DO
  END DO

  IF(mp_opt > 0) THEN
    CALL mpsendrecv2dew(tem9,nx,ny,nz,1,1,0,tem1)
    CALL mpsendrecv2dns(tem9,nx,ny,nz,1,1,0,tem1)
  END IF


  RETURN
END SUBROUTINE cal_div
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_divq                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Moist Conv.*1000. (g/kg/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_divq(tem9,u,v,qv,x,y,nx,ny,nz,tem1)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz), qv(nx,ny,nz)
  REAL :: x(nx), y(ny)
  REAL :: tem1(nx,ny,nz)  ! work array

  INTEGER :: i,j,k, istat

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=2,nz-1
    DO j=2,ny-1
      DO i=2,nx-1
        tem9(i,j,k)= -1. * 1000.0 * 1000.0 * 0.5 * (                    &
              ( u(i+1,j,k)*(qv(i,j,k)+qv(i+1,j,k))                      &
              -u(i,j,k)*(qv(i-1,j,k)+qv(i,j,k)) ) /(x(i+1)-x(i))        &
              + ( v(i,j+1,k)*(qv(i,j,k)+qv(i,j+1,k))                    &
              -v(i,j,k)*(qv(i,j-1,k)+qv(i,j,k)) ) /(y(j+1)-y(j)) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      tem9(i,j,   1)=tem9(i,j,   2)
      tem9(i,j,nz-1)=tem9(i,j,nz-2)
    END DO
  END DO

  IF(mp_opt > 0) THEN
    CALL mpsendrecv2dew(tem9,nx,ny,nz,1,1,0,tem1)
    CALL mpsendrecv2dns(tem9,nx,ny,nz,1,1,0,tem1)
  END IF

  RETURN
END SUBROUTINE cal_divq

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vtp                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate perturbation wind vectors.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vtp(tem7,tem8,tem9,uprt,vprt,wprt,nx,ny,nz,              &
           vtpunits,label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz), tem9(nx,ny,nz)
  REAL :: uprt(nx,ny,nz),vprt(nx,ny,nz),wprt(nx,ny,nz)
  INTEGER :: vtpunits,length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j,k, onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  onvf = 0
  CALL avgx(uprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  CALL avgy(vprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  CALL avgz(wprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)
  IF(vtpunits == 1) THEN
    label = '(m/s)'
    length =5
  ELSE IF(vtpunits == 2) THEN
    label = '(kts)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*1.943844
          tem8(i,j,k) = tem8(i,j,k)*1.943844
          tem9(i,j,k) = tem9(i,j,k)*1.943844
        END DO
      END DO
    END DO
  ELSE IF(vtpunits == 3) THEN
    label = '(MPH)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*2.236936
          tem8(i,j,k) = tem8(i,j,k)*2.236936
          tem9(i,j,k) = tem9(i,j,k)*2.236936
        END DO
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE cal_vtp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vtrstrm              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate wind streamline
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vtrstrm(tem7,tem8,tem9,u,v,w,nx,ny,nz,aspratio)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz),tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  REAL :: aspratio

  INTEGER :: i,j,k,onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  onvf = 0
  CALL avgx(u , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  CALL avgy(v , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  CALL avgz(w , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k)=aspratio*tem9(i,j,k)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_vtrstrm
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vtpstrm              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate the perturbation of wind streamlins.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vtpstrm(tem7,tem8,tem9,uprt,vprt,wprt,nx,ny,nz,          &
           aspratio)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz),tem9(nx,ny,nz)
  REAL :: uprt(nx,ny,nz),vprt(nx,ny,nz),wprt(nx,ny,nz)
  REAL :: aspratio

  INTEGER :: i,j,k,onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  onvf = 0
  CALL avgx(uprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  CALL avgy(vprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  CALL avgz(wprt , onvf,                                                &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)


  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k)=aspratio*tem9(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_vtpstrm

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_ vs                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vertical wind shear*1000(1/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vs(tem9,u,v,zp,tem7,tem8,nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: zp(nx,ny,nz),u(nx,ny,nz), v(nx,ny,nz)
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz)

  INTEGER :: i,j,k, onvf
  REAL :: tmp1, tmp2, tmp3
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  onvf = 0
  CALL avgx(u , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  onvf = 0
  CALL avgy(v , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tmp1 = ( zp(i,j,k+2) + zp(i,j,k+1) -                            &
                       zp(i,j,k) - zp(i,j,k-1) )*0.5
        tmp2 = tem7(i,j,k+1) - tem7(i,j,k-1)
        tmp3 = tem8(i,j,k+1) - tem8(i,j,k-1)
        tem9(i,j,k) = 1000.*SQRT((tmp2/tmp1)**2+(tmp3/tmp1)**2)

      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_vs

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_gric                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Richardson Number
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_gric(tem9,u,v,zp,pt,tem7,tem8,nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz)
  REAL :: zp(nx,ny,nz),pt(nx,ny,nz),u(nx,ny,nz), v(nx,ny,nz)
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz)

  INTEGER :: i,j,k, onvf
  REAL :: tmp1, tmp2, tmp3, tmp4
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  onvf = 0
  CALL avgx(u , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  onvf = 0
  CALL avgy(v , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tmp1 = ( zp(i,j,k+2) + zp(i,j,k+1) -                            &
                       zp(i,j,k) - zp(i,j,k-1) )*0.5
        tmp2 = tem7(i,j,k+1) - tem7(i,j,k-1)
        tmp3 = tem8(i,j,k+1) - tem8(i,j,k-1)
        tmp4 =  (tmp2/tmp1)**2 + (tmp3/tmp1)**2
        tmp2 = pt(i,j,k+1)-pt(i,j,k-1)
        tmp3 = 9.8*tmp2/pt(i,j,k)/tmp1/(tmp4+1.e-20)
        tem9(i,j,k) = SIGN( MIN(ABS(tmp3),10.0), tmp3 )
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_gric

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_avor                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate absolute Vort*10^5 (1/s)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!  Yunheng Wang(12/20/2002)
!  Added code for message passing version.
!
!-----------------------------------------------------------------------

SUBROUTINE cal_avor(tem9,u,v,x,y,nx,ny,nz,mode,flagsin,omega,           &
           sinlat,tem1,tem2,tem3)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz),sinlat(nx,ny)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)
  REAL :: tem1(nx,ny,nz), tem2(nx,ny,nz), tem3(nx,ny,nz)
  REAL :: omega
  INTEGER :: mode,flagsin

  INTEGER :: i,j,k
  REAL :: tmp1

  INCLUDE 'mp.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        tem9(i,j,k)= 1.0E5*(                                            &
               (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
               (4*(x(i+1)-x(i)))-                                       &
               (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
               (4*(y(j+1)-y(j))) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      tem9(i,j,   1)=tem9(i,j,   2)
      tem9(i,j,nz-1)=tem9(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      tem9(   1,j,k)=tem9(   2,j,k)
      tem9(nx-1,j,k)=tem9(nx-2,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      tem9(i,   1,k)=tem9(i,   2,k)
      tem9(i,ny-1,k)=tem9(i,ny-2,k)
    END DO
  END DO

  IF(mode == 1.OR.mode == 4.OR.mode == 6.OR.mode == 7) THEN
    IF( flagsin == 0) THEN
      CALL gtsinlat(nx,ny,x,y,sinlat,tem1,tem2, tem3)
      tmp1 = 2.0* omega
      DO j=1,ny
        DO i=1,nx
          sinlat(i,j) = tmp1 * sinlat(i,j)*1.0E5
        END DO
      END DO
      flagsin=1
    END IF
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem9(i,j,k) = tem9(i,j,k) + sinlat(i,j)
        END DO
      END DO
    END DO
  END IF

  IF(mp_opt > 0) THEN
    CALL mpsendrecv2dew(tem9,nx,ny,nz,1,1,0,tem1)
    CALL mpsendrecv2dns(tem9,nx,ny,nz,1,1,0,tem2)
  END IF

  RETURN
END SUBROUTINE cal_avor

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_qt                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Total water & vapor (g/kg)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_qt(tem9,qv,qscalar,nx,ny,nz)

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
  REAL    :: tem9(nx,ny,nz)
  REAL    :: qv(nx,ny,nz)
  REAL    :: qscalar(nx,ny,nz,nscalar)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem9(i,j,k) = qv(i,j,k)
        IF (P_QC > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QC)
        IF (P_QR > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QR)
        IF (P_QI > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QI)
        IF (P_QS > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QS)
        IF (P_QG > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QG)
        IF (P_QH > 0) tem9(i,j,k) = tem9(i,j,k) + qscalar(i,j,k,P_QH)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_qt


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_rhi                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate rhi value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_rhi(tem9,pt, pprt ,pbar,qv,tem1,tem2,nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem9(nx,ny,nz), tem1(nx,ny,nz), tem2(nx,ny,nz)
  REAL :: pt(nx,ny,nz), pprt(nx,ny,nz), pbar(nx,ny,nz), qv(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL temper (nx,ny,nz,pt, pprt ,pbar,tem1)
  CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, pbar,tem1,tem2)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem9(i,j,k) = qv(i,j,k)/tem2(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_rhi

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_xuv                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate wind vector (xuv)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_xuv(tem7,tem8,tem9,u,v,w,nx,ny,nz,xuvunits,              &
           label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz),tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  INTEGER :: xuvunits,length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j,k,onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  onvf = 0
  CALL avgx(u , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  CALL avgy(v , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  CALL avgz(w , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)

  IF(xuvunits == 1) THEN
    label = '(m/s)'
    length =5
  ELSE IF(xuvunits == 2) THEN
    label = '(kts)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*1.943844
          tem8(i,j,k) = tem8(i,j,k)*1.943844
          tem9(i,j,k) = tem9(i,j,k)*1.943844
        END DO
      END DO
    END DO
  ELSE IF(xuvunits == 3) THEN
    label = '(MPH)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*2.236936
          tem8(i,j,k) = tem8(i,j,k)*2.236936
          tem9(i,j,k) = tem9(i,j,k)*2.236936
        END DO
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE cal_xuv

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vtr                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate wind vector
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vtr(tem7,tem8,tem9,u,v,w,nx,ny,nz,vtrunits,label,        &
           length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz),tem9(nx,ny,nz)
  REAL :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  INTEGER :: vtrunits,length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j,k,onvf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  onvf = 0
  CALL avgx(u , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem7)
  CALL avgy(v , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem8)
  CALL avgz(w , onvf,                                                   &
           nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)

  IF(vtrunits == 1) THEN
    label = '(m/s)'
    length =5
  ELSE IF(vtrunits == 2) THEN
    label = '(kts)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*1.943844
          tem8(i,j,k) = tem8(i,j,k)*1.943844
          tem9(i,j,k) = tem9(i,j,k)*1.943844
        END DO
      END DO
    END DO
  ELSE IF(vtrunits == 3) THEN
    label = '(MPH)'
    length = 5
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem7(i,j,k) = tem7(i,j,k)*2.236936
          tem8(i,j,k) = tem8(i,j,k)*2.236936
          tem9(i,j,k) = tem9(i,j,k)*2.236936
        END DO
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE cal_vtr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vtrobs               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate wind observation value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vtrobs(dd,ff,drot, vtrunits)

  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp, obs_marksz

  INTEGER :: nobs
  COMMON /sfc_obs1/ nobs
  REAL :: latob(mxsfcob),lonob(mxsfcob)
  REAL :: obs1(mxsfcob),obs2(mxsfcob)
  COMMON /sfc_obs2/ latob,lonob,obs1,obs2

  REAL :: dd(mxsfcob),ff(mxsfcob)
  REAL :: drot

  INTEGER :: iob, vtrunits
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO iob=1,nobs
    IF(dd(iob) >= 0. .AND. dd(iob) < 360. .AND.                         &
          ff(iob) >= 0. .AND. ff(iob) < 98.) THEN
      CALL ddrotuv(1,lonob(iob),dd(iob),ff(iob),                        &
                    drot,obs1(iob),obs2(iob))
      obs1(iob)=0.51444*obs1(iob)
      obs2(iob)=0.51444*obs2(iob)
    ELSE
      obs1(iob)=-999.
      obs2(iob)=-999.
    END IF
  END DO
  obsset=1

  IF(vtrunits == 2) THEN    !! kts
    DO iob=1,nobs
      IF(obs1(iob) /= -999.) obs1(iob)= obs1(iob)*1.943844
      IF(obs2(iob) /= -999.) obs2(iob)= obs2(iob)*1.943844
    END DO
  ELSE IF(vtrunits == 3) THEN   !! MPH
    DO iob=1,nobs
      IF(obs1(iob) /= -999.) obs1(iob) = obs1(iob)*2.236936
      IF(obs1(iob) /= -999.) obs2(iob) = obs2(iob)*2.236936
    END DO
  END IF
  RETURN
END SUBROUTINE cal_vtrobs

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_viqc                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate viqc
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_viqc(tem7, qc, rhobar, zp, nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qc(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)+qc(i,j,k)*rhobar(i,j,k)               &
                        *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_viqc


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_viqr                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate viqr
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_viqr(tem7, qr, rhobar, zp, nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qr(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)+qr(i,j,k)*rhobar(i,j,k)               &
                        *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_viqr


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_viqi                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate viqi
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_viqi(tem7, qi, rhobar, zp, nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qi(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)+qi(i,j,k)*rhobar(i,j,k)               &
                        *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_viqi

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_viqs                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate viqs
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_viqs(tem7, qs, rhobar, zp, nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qs(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)+qs(i,j,k)*rhobar(i,j,k)               &
                        *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_viqs

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_viqh                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate viqh
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_viqh(tem7, qh, rhobar, zp, nx,ny,nz)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qh(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)+qh(i,j,k)*rhobar(i,j,k)               &
                        *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_viqh

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vil                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vert. Integ Liquid (kg/m2)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vil(tem7,qc,qr,rhobar,zp, nx,ny,nz,tem6)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qc(nx,ny,nz),qr(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem6(i,j,k) = qc(i,j,k) + qr(i,j,k)
        tem7(i,j,1) = tem7(i,j,1)+tem6(i,j,k)*rhobar(i,j,k)             &
                  *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_vil


!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vii                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vert. Integrated ice (kg/m2)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vii(tem7,qi,qs,qh,rhobar,zp, nx,ny,nz,tem6)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qi(nx,ny,nz), qs(nx,ny,nz), qh(nx,ny,nz)
  REAL :: rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem6(i,j,k) = qi(i,j,k) + qs(i,j,k) + qh(i,j,k)
        tem7(i,j,1) = tem7(i,j,1)+tem6(i,j,k)*rhobar(i,j,k)             &
                  *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE cal_vii

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vic                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate vic
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vic(tem2,qscalar,rhobar,zp,nx,ny,nz,tem1)

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
  REAL    :: tem2(nx,ny)
  REAL    :: qscalar(nx,ny,nz,nscalar)
  REAL    :: rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL    :: tem1(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem2(i,j)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k) = 0.0
        IF (P_QC > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QC)
        IF (P_QR > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QR)
        IF (P_QI > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QI)
        IF (P_QS > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QS)
        IF (P_QG > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QG)
        IF (P_QH > 0) tem1(i,j,k) = tem1(i,j,k) + qscalar(i,j,k,P_QH)
        tem2(i,j)   = tem2(i,j)+tem1(i,j,k)*rhobar(i,j,k)               &
                                           *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_vic

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_vit                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate vit ( vertically intergrated total water(kg/m**2))
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  5/10/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_vit(tem7,qv,qscalar,rhobar,zp,nx,ny,nz,tem6)

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: nx,ny,nz
  REAL    :: tem7(nx,ny)
  REAL    :: qscalar(nx,ny,nz,nscalar)
  REAL    :: qv(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL    :: tem6(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem6(i,j,k) = qv(i,j,k)
        IF (P_QC > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QC)
        IF (P_QR > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QR)
        IF (P_QI > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QI)
        IF (P_QS > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QS)
        IF (P_QG > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QG)
        IF (P_QH > 0) tem6(i,j,k) = tem6(i,j,k) + qscalar(i,j,k,P_QH)

        tem7(i,j)   = tem7(i,j)+tem6(i,j,k)*rhobar(i,j,k)               &
                                           *(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cal_vit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_pw                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cal_pw(tem7,qv,rhobar,zp,nx,ny,nz,tem6)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate pw ( precipitable water vapor(cm)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  5/10/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: qv(nx,ny,nz), rhobar(nx,ny,nz), zp(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx
        tem6(i,j,k) = qv(i,j,k)
        tem7(i,j,1) = tem7(i,j,1)                                       &
                    + tem6(i,j,k)*rhobar(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

! change kg/m**2 to cm

  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1) = 0.1*tem7(i,j,1)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_pw
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_tpr                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate tpr ( total precipatation rate )
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  5/15/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_tpr(tem7,prcrate,nx,ny,nz,tprunits,label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: prcrate(nx,ny)
  INTEGER :: tprunits, length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  IF(tprunits == 1) THEN
    label = 'Total precip. rate(mm/h)'
    length =24
    DO j=1,ny
      DO i=1,nx
           !RLC 1998/05/18
        !tem7(i,j,1) = prcrate(i,j)/100./3600. !!kg/(m**2*s) -> mm/h
        tem7(i,j,1) = prcrate(i,j)*3600.   !!kg/(m**2*s) -> mm/h
      END DO
    END DO
  ELSE IF(tprunits == 2) THEN   !! kg/(m**2*s) -> in/h
    label = 'Total precip. rate(in/h)'
    length =24
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = prcrate(i,j)*0.039370079*3600.
      END DO
    END DO
  ELSE
    WRITE(*,'(1x,a,I0)') 'ERROR: Wrong namelist parameter tprunits = ',tprunits
    CALL arpsstop('ERROR: Wrong namelist parameter tprunits.',1)
  END IF

  RETURN
END SUBROUTINE cal_tpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_gpr                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate gpr ( Grid-scale precip. rate)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  5/15/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_gpr(tem7,prcrate,prcrate1,nx,ny,nz,                      &
           gprunits,label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: prcrate(nx,ny), prcrate1(nx,ny)
  INTEGER :: gprunits, length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  IF(gprunits == 1) THEN
    label = 'Grid-scale precip. rate(mm/h)'
    length =29
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = (prcrate(i,j)+prcrate1(i,j))*3600. !kg/(m**2*s) -> mm/h
      END DO
    END DO
  ELSE IF(gprunits == 2) THEN
                          !! kg/(m**2*s) -> in/h
    label = 'Grid-scale precip. rate(in/h)'
    length =29
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = (prcrate(i,j)+prcrate1(i,j)) *0.039370079*3600.
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE cal_gpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_cpr                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate cpr ( Convective precip. rate )
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  5/15/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_cpr(tem7,prcrate,nx,ny,nz,cprunits,label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz)
  REAL :: prcrate(nx,ny)
  INTEGER :: cprunits, length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx
      tem7(i,j,1)=0.
    END DO
  END DO

  IF(cprunits == 1) THEN
    label = 'Convective precip. rate(mm/h)'
    length =29
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = prcrate(i,j)*3600. !!kg/(m**2*s) -> mm/h
      END DO
    END DO
  ELSE IF(cprunits == 2) THEN
                          !! kg/(m**2*s) -> in/h
    label = 'Convective precip. rate(in/h)'
    length =29
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = prcrate(i,j)*0.039370079*3600.
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE cal_cpr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE cal_strm                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate strom motion
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_strm(tem7, tem8,ustrm,vstrm,strmunits,nx,ny,nz,          &
           label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: tem7(nx,ny,nz), tem8(nx,ny,nz)
  REAL :: ustrm(nx,ny) ,vstrm(nx,ny)
  INTEGER :: strmunits, length
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny
    DO i=1,nx-1
      tem7(i,j,1)=(ustrm(i,j)+ustrm(i+1,j))*0.5
    END DO
  END DO
  DO j=1,ny-1
    DO i=1,nx
      tem8(i,j,1)=(vstrm(i,j)+vstrm(i,j+1))*0.5
    END DO
  END DO

  IF(strmunits == 1) THEN
    label = 'storm motion (m/s)'
    length =18
  ELSE IF(strmunits == 2) THEN
    label = 'storm motion (kts)'
    length =18
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)*1.943844
        tem8(i,j,1) = tem8(i,j,1)*1.943844
      END DO
    END DO
  ELSE IF(strmunits == 3) THEN
    label = 'storm motion (MPH)'
    length = 18
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = tem7(i,j,1)*2.236936
        tem8(i,j,1) = tem8(i,j,1)*2.236936
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE cal_strm
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE cal_arbvtr                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate arbitrary vector
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  6/1/2007
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_arbvtr(tem7,tem8,arbu,arbv,vtraunits,nx,ny,nz,          &
           varname,label,length)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL    :: tem7(nx,ny,nz), tem8(nx,ny,nz)
  REAL    :: arbu(nx,ny) ,arbv(nx,ny)
  INTEGER :: vtraunits, length
  CHARACTER (LEN=6) :: varname
  CHARACTER (LEN=*) :: label

  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF(vtraunits == 1) THEN
    write(label,'(a6,a6)') varname,' (m/s)'
    length =18
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = arbu(i,j)
        tem8(i,j,1) = arbv(i,j)
      END DO
    END DO
  ELSE IF(vtraunits == 2) THEN
    write(label,'(a6,a6)') varname,' (kts)'
    length =12
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = arbu(i,j)*1.943844
        tem8(i,j,1) = arbv(i,j)*1.943844
      END DO
    END DO
  ELSE IF(vtraunits == 3) THEN
    write(label,'(a6,a6)') varname,' (MPH)'
    length = 12
    DO j=1,ny
      DO i=1,nx
        tem7(i,j,1) = arbu(i,j)*2.236936
        tem8(i,j,1) = arbv(i,j)*2.236936
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE cal_arbvtr

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_dist                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_dist(haxisu,dx,dy,x01,y01,x02,y02,slicopt,                  &
           tmpx,tmpy,distc)


  IMPLICIT NONE
  INTEGER :: haxisu, slicopt
  REAL :: dx,dy
  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  REAL :: x01,y01                  ! the first  point of interpolation
  REAL :: x02,y02                  ! the second point of interpolation
  REAL :: tmpx, tmpy
  CHARACTER (LEN=*) :: distc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF ( slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR.   &
       slicopt == 10 .OR. slicopt == 11 .OR. slicopt==12 ) THEN
    IF (haxisu == 0) THEN
      tmpx = dx          !!!km
      tmpy = dy         !!!km
      x101 = x01
      y101 = y01
      x102 = x02
      y102 = y02
      WRITE(distc,'('' KM'')')
    ELSE IF(haxisu == 1 ) THEN
      tmpx = dx*0.62137 !!! mile
      tmpy = dy*0.62137 !!! mile
      x101 = x01*0.62137
      y101 = y01*0.62137
      x102 = x02*0.62137
      y102 = y02*0.62137
      WRITE(distc,'('' MILE'')')
    ELSE IF(haxisu == 2) THEN
      tmpx = dx*0.53997 !!!naut mile
      tmpy = dy*0.53997 !!!naut mile
      x101 = x01*0.53997
      y101 = y01*0.53997
      x102 = x02*0.53997
      y102 = y02*0.53997
      WRITE(distc,'('' NAUT MILE'')')
    ELSE IF(haxisu == 3) THEN
      tmpx = dx*3.28084 !!!kft
      tmpy = dy*3.28084 !!!kft
      x101 = x01*3.28084
      y101 = y01*3.28084
      x102 = x02*3.28084
      y102 = y02*3.28084
      WRITE(distc,'('' KFT'')')
    END IF
  END IF
!  IF(mode.eq.5) THEN
!    length=len(distc)
!    CALL strmin(distc,length)
!    write(distc,'(a)') '('//distc(2:length)
!  ENDIF
  RETURN
END SUBROUTINE cal_dist

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE cal_Nx                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cal_Nx(nx,ny,nz,rhobar,tz,qscalar,temscalar)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Diagnose Number concentration of each microphysics species
!  given the single-moment microphysics scheme using the appropriate
!  intercept parameter, mixing ratio, and densities of each species
!  Currently works for the ARPS Lin scheme, as well as the fixed
!  and Diagnostic N0r versions of the WSM6 scheme and the single-moment
!  version of the Milbrandt and Yau (2005) multi-moment scheme
!
!  For now the calculation assumes that the size distributions are
!  exponential (i.e. alpha (shape parameter) = 0 in the gamma distribution
!  formula).
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  2/26/07
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'

  INTEGER :: nx,ny,nz
  REAL :: rhobar(nx,ny,nz)
  REAL :: tz(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: temscalar(nx,ny,nz,nscalar)

  REAL,PARAMETER :: pi = 3.14159265

  REAL :: N0r, N0s, N0g, N0h
  REAL :: rhor,rhoi,rhos,rhog,rhoh

  INTEGER :: i,j,k,nq

  rhor = 1000.

  IF(mphyopt == 2) THEN    ! Lin scheme
    N0r = n0rain
    N0s = n0snow
    N0h = n0hail
    rhos = rhosnow
    rhoh = rhohail
!    rhos = 100.
!    rhoh = 913.

  ELSE IF(mphyopt == 5) THEN    ! original WSM6 scheme
    N0r = 8.0e6
    N0g = 4.0e6
    rhog = 500.
    rhos = 100.
  ELSE IF(mphyopt == 7) THEN    ! diagnostic N0r WSM6 scheme
    N0g = 4.0e6
    rhog = 500.
    rhos = 100.
  ELSE IF(mphyopt == 8) THEN
    N0r = n0rain
    N0s = n0snow
    N0g = n0grpl
    N0h = n0hail

    rhoi = rhoice
    rhos = rhosnow
    rhog = rhogrpl
    rhoh = rhohail
  ELSE                     ! Schultz, Straka schemes : assumed same as default Lin for now
    N0r = 8.0e6
    N0s = 3.0e6
    N0h = 4.0e4
  END IF

  temscalar = 0.0

  DO nq=1,nscalar
   DO i=2,nx-1
    DO j=2,ny-1
      DO k=2,nz-1
        IF(qnames(nq) == 'qc') THEN
          IF(mphyopt == 2) THEN
            temscalar(i,j,k,nq) = 1.0e9   ! Constant
          ELSE IF(mphyopt == 5 .or. mphyopt == 7) THEN
            temscalar(i,j,k,nq) = 3.0e8   ! Constant
          ELSE IF(mphyopt == 8) THEN
            temscalar(i,j,k,nq) = 1.0e8   ! Constant
          END IF
        ELSE IF(qnames(nq) == 'qr') THEN
          IF(mphyopt == 7) THEN
            ! Diagnostic N0r
            N0r = 7835.5*1000*(qscalar(i,j,k,nq)*rhobar(i,j,k)*1000)**0.681
          END IF
          !print*,"Here inside cal_Nx for rain!"
          temscalar(i,j,k,nq) = (N0r**(3./4.))*(rhobar(i,j,k)*qscalar(i,j,k,nq)  &
                                /(pi*rhor))**(1./4.)
        ELSE IF(qnames(nq) == 'qi') THEN
          IF(mphyopt == 2) THEN
            temscalar(i,j,k,nq) = 0.0  ! Not sure how to calculate Ni for Lin scheme
          ELSE IF(mphyopt == 5 .or. mphyopt == 7) THEN  !Mass-based Ni for WSM6 scheme
            temscalar(i,j,k,nq) = 5.38e7*(rhobar(i,j,k)*qscalar(i,j,k,nq))**0.75
          ELSE IF(mphyopt == 8) THEN
            ! MY scheme appears to use the following Cooper equation for ice number
            ! concentration for the single-moment scheme.  Other tendencies are calculated
            ! but do not appear to be used in the single-moment case.
            temscalar(i,j,k,nq) = 5.*exp(0.304*(273.15-max(233.,tz(i,j,k))))
          END IF
        ELSE IF(qnames(nq) == 'qs') THEN
          IF(mphyopt == 5 .or. mphyopt == 7) THEN  ! Temperature-dependent N0s (and thus Ns)
            N0s = min(2e6*exp(0.12*(273.15-tz(i,j,k))),1e11)
          END IF
          temscalar(i,j,k,nq) = (N0s**(3./4.))*(rhobar(i,j,k)*qscalar(i,j,k,nq)  &
                                /(pi*rhos))**(1./4.)
        ELSE IF(qnames(nq) == 'qg') THEN
          temscalar(i,j,k,nq) = (N0g**(3./4.))*(rhobar(i,j,k)*qscalar(i,j,k,nq)  &
                                /(pi*rhog))**(1./4.)
        ELSE IF(qnames(nq) == 'qh') THEN
          temscalar(i,j,k,nq) = (N0h**(3./4.))*(rhobar(i,j,k)*qscalar(i,j,k,nq)  &
                                /(pi*rhoh))**(1./4.)
        END IF
      END DO
    END DO
   END DO
  END DO

END SUBROUTINE cal_Nx

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES setcords                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setcords(xl,xr,yb,yt,dx,dy, slicopt,                            &
           x1,x2,y1,y2,xlabel,ylabel,xstep,ystep,xmstep,ymstep)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Set coordinate related variables according to desired length unit
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!        xl,xr,yb,yt,dx,dy, slicopt
!  OUTPUT:
!        x1,x2,y1,y2,xlabel,ylabel,xstep,ystep,xmstep,ymstep
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  REAL,    INTENT(IN)  :: xl,xr,yb,yt,dx,dy
  INTEGER, INTENT(IN)  :: slicopt
  REAL,    INTENT(OUT) :: x1,x2,y1,y2
  REAL,    INTENT(OUT) :: xstep,ystep, xmstep, ymstep
  CHARACTER (LEN=*), INTENT(OUT) :: xlabel
  CHARACTER (LEN=*), INTENT(OUT) :: ylabel

  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF(slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.  &
     slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9) THEN

    IF(haxisu == 0) THEN
      x1=xl
      x2=xr
      y1=yb
      y2=yt
      !write(0,*) y1,y2,yb,yt   ! G95 strange
      WRITE(xlabel,'(a15)')'(km)'
      WRITE(ylabel,'(a15)')'(km)'
      !write(0,*) y1,y2,yb,yt   ! G95 changes the INTENT(IN) values
    ELSE IF(haxisu == 1 ) THEN
      x1=xl*0.62137
      x2=xr*0.62137
      y1=yb*0.62137
      y2=yt*0.62137
      WRITE(xlabel,'(a15)')'(mile)'
      WRITE(ylabel,'(a15)')'(mile)'
    ELSE IF(haxisu == 2) THEN
      x1=xl*0.53997
      x2=xr*0.53997
      y1=yb*0.53997
      y2=yt*0.53997
      WRITE(xlabel,'(a15)')'(naut mile)'
      WRITE(ylabel,'(a15)')'(naut mile)'
    ELSE IF(haxisu == 3) THEN
      x1=xl*3.28084
      x2=xr*3.28084
      y1=yb*3.28084
      y2=yt*3.28084
      WRITE(xlabel,'(a15)')'(kft)'
      WRITE(ylabel,'(a15)')'(kft)'
    END IF
    IF(tickopt == 0) THEN
      xstep = dx
      ystep = dy
      xmstep = 0.
      ymstep = 0.
    ELSE IF(tickopt == 1) THEN
      xstep = hmintick
      ystep = hmintick
      xmstep = hmajtick
      ymstep = hmajtick
    END IF

  END IF

  IF(slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR. slicopt==12) THEN
    IF(haxisu == 0) THEN
      x1=xl
      x2=xr
      WRITE(xlabel,'(a15)')'(km)'
    ELSE IF(haxisu == 1 ) THEN
      x1=xl*0.62137
      x2=xr*0.62137
      WRITE(xlabel,'(a15)')'(mile)'
    ELSE IF(haxisu == 2) THEN
      x1=xl*0.53997
      x2=xr*0.53997
      WRITE(xlabel,'(a15)')'(naut mile)'
    ELSE IF(haxisu == 3) THEN
      x1=xl*3.28084
      x2=xr*3.28084
      WRITE(xlabel,'(a15)')'(kft)'
    END IF
    IF(vaxisu == 0) THEN
      y1=yb
      y2=yt
      WRITE(ylabel,'(a15)')'(km)'
    ELSE IF(vaxisu == 1 ) THEN
      y1=yb*0.62137
      y2=yt*0.62137
      WRITE(ylabel,'(a15)')'(mile)'
    ELSE IF(vaxisu == 2) THEN
      y1=yb*0.53997
      y2=yt*0.53997
      WRITE(ylabel,'(a15)')'(naut mile)'
    ELSE IF(vaxisu == 3) THEN
      y1=yb*3.28084
      y2=yt*3.28084
      WRITE(ylabel,'(a15)')'(kft)'
    ELSE IF(vaxisu == 4) THEN
      y1=yb*3.28084
      y2=yt*3.28084
      WRITE(ylabel,'(a15)')'(presure)'
    END IF
    IF(tickopt == 0) THEN
      xstep = dx
      ystep = dy
      xmstep = 0.
      ymstep = 0.
    ELSE IF(tickopt == 1) THEN
      xstep = hmintick
      ystep = vmintick
      xmstep = hmajtick
      ymstep = vmajtick
    END IF
  END IF

  IF(slicopt == 10 .OR. slicopt == 11) THEN
    IF(haxisu == 0) THEN
      x1=xl
      x2=xr
      WRITE(xlabel,'(a15)')'(km)'
    ELSE IF(haxisu == 1 ) THEN
      x1=xl*0.62137
      x2=xr*0.62137
      WRITE(xlabel,'(a15)')'(mile)'
    ELSE IF(haxisu == 2) THEN
      x1=xl*0.53997
      x2=xr*0.53997
      WRITE(xlabel,'(a15)')'(naut mile)'
    ELSE IF(haxisu == 3) THEN
      x1=xl*3.28084
      x2=xr*3.28084
      WRITE(xlabel,'(a15)')'(kft)'
    END IF
    IF(vaxisu == 0) THEN
      y1=yb
      y2=yt
      WRITE(ylabel,'(a15)')'(cm)'
    END IF
    IF(tickopt == 0) THEN
      xstep = dx
      ystep = dy
      xmstep = 0.
      ymstep = 0.
    ELSE IF(tickopt == 1) THEN
      xstep = hmintick
      ystep = vmintick
      xmstep = hmajtick
      ymstep = vmajtick
    END IF
  END IF

  RETURN
END SUBROUTINE setcords
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_CONTOUR                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_contour(ncont, tcont)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set-up ccontour values for when plot option set to 11 .
!    right now only work for several variables:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!
!  MODIFICATION HISTORY:
!  2/18/97
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!  tcon     the total number of contours.
!  ncon     the contour array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arpsplt.inc'

  INCLUDE 'mp.inc'

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  INTEGER :: setcontopt, setcontnum
  CHARACTER (LEN=12) :: setcontvar(maxuneva)
  REAL :: setconts(maxunevm,maxuneva)
  COMMON /setcont_var/setcontvar
  COMMON /setcon_par/setcontopt,setcontnum,setconts

  INTEGER :: i, LEN, var, ncont
  REAL    :: tcont(maxunevm)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ncont = 0
  var = 0

  LEN=12
  CALL strlnth( varname, LEN)
  IF(setcontopt == 1) THEN
    DO i = 1,setcontnum
      IF(varname(1:LEN) == setcontvar(i)(1:LEN)) THEN
        var = i
        GO TO 10
      END IF
    END DO
  END IF

  IF(var == 0) RETURN

  10    CONTINUE

  DO i=1,maxunevm
    IF( setconts(i,var) == -9999.0) THEN
      GO TO 100
    END IF
  END DO
  100   ncont = i-1

  DO i=1,ncont
    tcont(i) = setconts(i,var)
  END DO

  IF(var /= 0 .AND. myproc == 0)     &
    WRITE(6,'(1x,a,I3,a,32F12.5)') '* irregular contours - ',           &
                               ncont, ':',(tcont(i),i=1,ncont)

  RETURN
END SUBROUTINE get_contour
