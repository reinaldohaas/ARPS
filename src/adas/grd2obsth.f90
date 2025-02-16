










!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRDTOSNG                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grdtosng(nx,ny,nxlg,nylg,nz,nz_tab,mxsng,nvar,nobsng,        &
           xs,ys,xslg,yslg,zp,icatg,anx,qback,hqback,nlvqback,          &
           shght,su,sv,spres,stheta,sqv,                                &
           stnsng,isrcsng,icatsng,elevsng,xsng,ysng,                    &
           obsng,qobsng,qualsng,odifsng,oanxsng,thesng,trnsng,          &
           indexsng)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Gridded data are interpolated to single-level locations to
!  produce observation differences.  For now, bilinear interpolation
!  is used to create a vertical sounding at the obs location
!  and that is interpolated in the vertical to obtain the
!  bacground value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, July, 1995
!
!  MODIFICATION HISTORY:
!  Sept, 1996 (KB)
!  Added count-reporting diagnostic output.
!
!  Dec., 1998 (KB)
!  Added correlation categories.
!
!  Feb., 1999 (KB)
!  Moved calculation of oanxsng variables to outside the IF (good)
!  block, so that it will always be available, for future
!  calculations.   Needed for determining qv thresholding from RH
!  thresholding.   Added sychronization of wind error flags, if one
!  component is judged to be in error, the other is also.
!
!  September, 2002 (KB)
!  Added processing for trnsng
!
!  December, 2005 (KWT)
!  MPI version.
!
!  November, 2005 (KWT)
!  Add MPI for radial velocity.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'

!
  INTEGER :: nx,ny,nxlg,nylg,nz,nz_tab
  INTEGER :: mxsng,nvar
  INTEGER :: nobsng
!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx), xslg(nxlg)
  REAL :: ys(ny), yslg(nylg)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
  INTEGER :: icatg(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Background fields
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
  REAL :: qback(nvar,nz_tab)
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback
!
!-----------------------------------------------------------------------
!
!  Column-interpolated values (used internally only)
!
!-----------------------------------------------------------------------
!
  REAL :: shght(nz)
  REAL :: su(nz)
  REAL :: sv(nz)
  REAL :: spres(nz)
  REAL :: stheta(nz)
  REAL :: sqv(nz)
!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnsng(mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)

  REAL :: elevsng(mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: obsng(nvar,mxsng)
  REAL :: qobsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  REAL :: odifsng(nvar,mxsng)
  REAL :: oanxsng(nvar,mxsng)
  REAL :: thesng(mxsng)
  REAL :: trnsng(mxsng)
  INTEGER :: indexsng(mxsng)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,kgrid,ktab,kntdom
  INTEGER :: ipt,jpt,ireturn,nlevs,k
  REAL :: selev,hgtgnd,hgtagl
  REAL :: whigh,wlow,wqhigh,wqlow,qbthe
!
!-----------------------------------------------------------------------
!
!  Beginning of executable code
!
!-----------------------------------------------------------------------
!
  kntdom=0
  DO ista=1,nobsng
!
!-----------------------------------------------------------------------
!
!   Combined ("suprob") obs are already marked to never be used.
!
!-----------------------------------------------------------------------
!
    IF (isrcsng(ista) == -1) CYCLE

    trnsng(ista) = -999.0
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL findlc(nx,ny,xs,ys,xsng(ista),ysng(ista),                      &
                ipt,jpt,ireturn)
!
!-----------------------------------------------------------------------
!
!  Interpolated gridded data if extrapolation is not indicated.
!
!  The extrapolation code isn't active.  DON'T CONSIDER CHANGING THIS!
!  It will break MPI!!!
!
!-----------------------------------------------------------------------
!
    IF(ireturn >= 0) THEN
!      write(6,'(a,a)')  stnsng(ista),
!    :            ' inside ARPS domain, processing'

!
!-----------------------------------------------------------------------
!
!  In MPI mode, some obs will fall on more than once processor.  Make
!  sure no ob gets counted twice.
!
!-----------------------------------------------------------------------
!
!      IF (mp_opt == 0 .OR. indexsng(ista) == myproc)                    &
      IF ( indexsng(ista) == myproc) kntdom=kntdom+1

!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6)')
!    :  ' location x= ',(0.001*xsng(ista)),
!    :         '   y= ',(0.001*ysng(ista)),
!    :  ' found near pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2,/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
      CALL colint(nx,ny,nz,nvar,                                        &
             xs,ys,zp,xsng(ista),ysng(ista),ipt,jpt,anx,                &
             su,sv,stheta,spres,shght,sqv,selev,                        &
             nlevs)
!
!-----------------------------------------------------------------------
!
!  Assign category
!
!-----------------------------------------------------------------------
!
      icatsng(ista)=icatg(ipt,jpt)
!
!-----------------------------------------------------------------------
!
!  Interpolate in the vertical at each obs level and compute
!  obs differences.
!
!-----------------------------------------------------------------------
!
      hgtgnd=0.5*(shght(2)+shght(1))
      trnsng(ista)=hgtgnd
      DO kgrid=2,nlevs-2
        IF(shght(kgrid) > elevsng(ista)) EXIT
      END DO
      whigh=(elevsng(ista)-shght(kgrid-1))/                             &
             (shght(kgrid)-shght(kgrid-1))
      wlow=1.-whigh
!
      hgtagl=AMIN1(hqback(nlvqback),                                    &
                   AMAX1(hqback(1),(elevsng(ista)-hgtgnd)))
      DO ktab=2,nlvqback-1
        IF(hqback(ktab) > hgtagl) EXIT
      END DO
      wqhigh=(hgtagl-hqback(ktab-1))/                                   &
             (hqback(ktab)-hqback(ktab-1))
      wqlow=1.-wqhigh
!
      oanxsng(1,ista)=(wlow*su(kgrid-1)+whigh*su(kgrid))
      IF(qualsng(1,ista) > 0) THEN
        odifsng(1,ista)=obsng(1,ista)-oanxsng(1,ista)
        qobsng(1,ista)=qobsng(1,ista)/                                  &
              (wqlow*qback(1,ktab-1)+wqhigh*qback(1,ktab))
      END IF

      oanxsng(2,ista)=(wlow*sv(kgrid-1)+whigh*sv(kgrid))
      IF(qualsng(2,ista) > 0) THEN
        odifsng(2,ista)=obsng(2,ista)-oanxsng(2,ista)
        qobsng(2,ista)=qobsng(2,ista)/                                  &
              (wqlow*qback(2,ktab-1)+wqhigh*qback(2,ktab))
      END IF

      oanxsng(3,ista)=EXP( wlow*LOG(spres(kgrid-1))+                    &
                            whigh*LOG(spres(kgrid)) )
      IF(qualsng(3,ista) > 0) THEN
        odifsng(3,ista)=obsng(3,ista)-oanxsng(3,ista)
        qobsng(3,ista)=qobsng(3,ista)/                                  &
              (wqlow*qback(3,ktab-1)+wqhigh*qback(3,ktab))
      END IF

      oanxsng(4,ista)=(wlow*stheta(kgrid-1)+whigh*stheta(kgrid))
      IF(qualsng(4,ista) > 0) THEN
        odifsng(4,ista)=obsng(4,ista)-oanxsng(4,ista)
        qbthe=wqlow*qback(4,ktab-1)+wqhigh*qback(4,ktab)
        thesng(ista)=( (obsng(4,ista)/qobsng(4,ista))+                  &
                       (oanxsng(4,ista)/qbthe) ) /                      &
                       ((1./qobsng(4,ista)) + (1./qbthe))
!        print *, ' Theta:anx,obs,qsrc,qback: ',
!    :            oanxsng(4,ista),obsng(4,ista),qobsng(4,ista),qbthe
        qobsng(4,ista)=qobsng(4,ista)/qbthe
!        print *, ' theta-sng qobs normalized: ',
!    :            thesng(ista),qobsng(4,ista)
      ELSE
        thesng(ista)=(wlow*stheta(kgrid-1)+whigh*stheta(kgrid))
      END IF

      oanxsng(5,ista)=(wlow*sqv(kgrid-1)+whigh*sqv(kgrid))
      IF(qualsng(5,ista) > 0) THEN
        odifsng(5,ista)=obsng(5,ista)-oanxsng(5,ista)
        qobsng(5,ista)=qobsng(5,ista)/                                  &
              (wqlow*qback(5,ktab-1)+wqhigh*qback(5,ktab))
      END IF
!
!
    ELSE
!
!      write(6,'(a,a)')  stnsng(ista),
!    :            ' outside ARPS domain, skipping'
!

!
!   The rule for "is an observation to be used" is if the ob is in the domain
!   somewhere.  For MPI, an ob that we might need might be outside of the local
!   processor domain, but inside of the large domain.  Make sure it stays
!   valid.
!

!      ireturn = -1
!      IF (mp_opt > 0 ) THEN
!          CALL findlc(nxlg,nylg,xslg,yslg,xsng(ista),ysng(ista),        &
!                ipt,jpt,ireturn)
!      END IF
!
!      IF ( ireturn < 0 ) THEN
       IF ( indexsng(ista) < 0) THEN   ! no one owns it
         isrcsng(ista)=0
         icatsng(ista)=1

!      write(6,'(a,a,i4)') ' Extrapolation status returned from',
!    :                        ' subroutine findlc: ',ireturn
!      IF(ireturn.eq.-1) write(6,'(a)')
!    :      ' Extrapolation in x detected'
!      IF(ireturn.eq.-2) write(6,'(a)')
!    :      ' Extrapolation in y detected'
!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6/)')
!    :  ' location x= ',(0.001*xsng(ista)),
!    :         '   y= ',(0.001*ysng(ista)),
!    :  ' nearest pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))

         DO ivar=1,nvar
           qualsng(ivar,ista)=-9
         END DO
       END IF
     END IF
  END DO

  IF (mp_opt == 0 ) THEN
    WRITE(6,'(//,1x,i5,a,i5,a)') kntdom,' of ',nobsng,                  &
          ' single-level stations are inside the ARPS domain'
  ELSE
    IF (myproc == 0) WRITE(6,'(//)')
      CALL flush(6)
      CALL mpbarrier
    DO k=0,nprocs-1
      IF (myproc == k .AND. kntdom > 0) THEN
        WRITE(6,'(1x,i5,a,i5,a,a,1x,i4)') kntdom,' of ',nobsng,         &
          ' single-level stations are inside the ARPS domain',          &
          ' for proc ',myproc
        CALL flush(6)
      END IF
      CALL mpbarrier
    END DO
  END IF

!  IF(.FALSE.) THEN    ! Do not do it by default. It may desired for 
!                        ! some case with 3dvar, then change it to .TRUR.
!    nobsng=kntdom
!    write(*,*) ' Change the observation number to real used data number',nobsng
!  END IF
  

  RETURN
END SUBROUTINE grdtosng
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRDTOUA                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE grdtoua(nx,ny,nxlg,nylg,nz,nz_tab,nzua,mxua,nvar,nobsua,     &
           xs,ys,xslg,yslg,zp,anx,qback,hqback,nlvqback,                &
           shght,su,sv,spres,stheta,sqv,                                &
           stnua,isrcua,elevua,xua,yua,hgtua,                           &
           obsua,qobsua,qualua,nlevsua,                                 &
           odifua,oanxua,theua,trnua,indexua)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Gridded data are interpolated to upper air locations to
!  produce observation differences.  For now, bilinear interpolation
!  is used to create a vertical sounding at the obs location
!  and that is interpolated in the vertical to obtain the
!  bacground value.
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, July, 1995
!
!  MODIFICATION HISTORY:
!  Sept, 1996 (KB)
!  Added count-reporting diagnostic output.
!
!  September, 2002 (KB)
!  Added processing for trnua
!
!  January, 2006 (KWT)
!  MPI update.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
!
  INTEGER :: nx,ny,nxlg,nylg,nz,nz_tab
  INTEGER :: nzua,mxua,nvar
  INTEGER :: nobsua
!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx), xslg(nxlg)
  REAL :: ys(ny), yslg(nylg)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
!-----------------------------------------------------------------------
!
!  Background fields
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
  REAL :: qback(nvar,nz_tab)
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback
!
  REAL :: shght(nz)
  REAL :: su(nz)
  REAL :: sv(nz)
  REAL :: spres(nz)
  REAL :: stheta(nz)
  REAL :: sqv(nz)
!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnua(mxua)
  INTEGER :: isrcua(mxua)
  REAL :: elevua(mxua)
  REAL :: xua(mxua),yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: theua(nzua,mxua)
  REAL :: trnua(mxua)
  INTEGER :: indexua(mxua)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,k,kgrid,ktab,kntdom
  INTEGER :: ipt,jpt,ireturn,nlevs
  REAL :: selev,hgtgnd,hgtagl
  REAL :: whigh,wlow,wqhigh,wqlow,qbthe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  kntdom=0
  DO ista=1,nobsua
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL findlc(nx,ny,xs,ys,xua(ista),yua(ista),                        &
                ipt,jpt,ireturn)
!
!-----------------------------------------------------------------------
!
!  Interpolate if extrapoltion is not indicated.
!
!-----------------------------------------------------------------------
!
    IF(ireturn >= 0) THEN
!      write(6,'(a,a)')  stnua(ista),
!    :            ' inside ARPS domain, processing'

!
!-----------------------------------------------------------------------
!
!  In MPI mode, some obs will fall on more than once processor.  Make
!  sure no ob gets counted twice.
!
!-----------------------------------------------------------------------
!
!      IF (mp_opt == 0 .OR. indexua(ista) == myproc)                    &
      IF (indexua(ista) == myproc) kntdom=kntdom+1

!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6)')
!    :  ' location x= ',(0.001*xua(ista)),
!    :         '   y= ',(0.001*yua(ista)),
!    :  ' found near pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2,/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
      CALL colint(nx,ny,nz,nvar,                                        &
             xs,ys,zp,xua(ista),yua(ista),ipt,jpt,anx,                  &
             su,sv,stheta,spres,shght,sqv,selev,                        &
             nlevs)
!
!-----------------------------------------------------------------------
!
!  Convert pressure to ln(p) for interpolation
!
!-----------------------------------------------------------------------
!
      DO k=1,nlevs
        spres(k)=LOG(spres(k))
      END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate in the vertical at each obs level and compute
!  obs differences.
!
!-----------------------------------------------------------------------
!
      hgtgnd=0.5*(shght(2)+shght(1))
      trnua(ista)=hgtgnd
      DO k=1,nlevsua(ista)
        DO kgrid=2,nlevs-1
          IF(shght(kgrid) > hgtua(k,ista)) GO TO 260
        END DO
        DO ivar=1,nvar
          qualua(ivar,k,ista)=-9
        END DO
        CYCLE
        260       CONTINUE
        whigh=(hgtua(k,ista)-shght(kgrid-1))/                           &
               (shght(kgrid)-shght(kgrid-1))
        wlow=1.-whigh
!
        hgtagl=AMIN1(hqback(nlvqback),                                  &
                     AMAX1(hqback(1),(hgtua(k,ista)-hgtgnd)))
        DO ktab=2,nlvqback-1
          IF(hqback(ktab) > hgtagl) EXIT
        END DO
!        281       CONTINUE
        wqhigh=(hgtagl-hqback(ktab-1))/                                 &
               (hqback(ktab)-hqback(ktab-1))
        wqlow=1.-wqhigh
!
        IF(qualua(1,k,ista) > 0) THEN
          oanxua(1,k,ista)=wlow*su(kgrid-1)+whigh*su(kgrid)
          odifua(1,k,ista)=obsua(1,k,ista)-oanxua(1,k,ista)
          qobsua(1,k,ista)=qobsua(1,k,ista)/                            &
              (wqlow*qback(1,ktab-1)+wqhigh*qback(1,ktab))
        END IF
        IF(qualua(2,k,ista) > 0) THEN
          oanxua(2,k,ista)=wlow*sv(kgrid-1)+whigh*sv(kgrid)
          odifua(2,k,ista)=obsua(2,k,ista)-oanxua(2,k,ista)
          qobsua(2,k,ista)=qobsua(2,k,ista)/                            &
              (wqlow*qback(2,ktab-1)+wqhigh*qback(2,ktab))
        END IF
        IF(qualua(3,k,ista) > 0) THEN
          oanxua(3,k,ista)=EXP( wlow*spres(kgrid-1)+                    &
                               whigh*spres(kgrid) )
          odifua(3,k,ista)=obsua(3,k,ista)-oanxua(3,k,ista)
          qobsua(3,k,ista)=qobsua(3,k,ista)/                            &
              (wqlow*qback(3,ktab-1)+wqhigh*qback(3,ktab))
        END IF
        IF(qualua(4,k,ista) > 0) THEN
          oanxua(4,k,ista)=wlow*stheta(kgrid-1)+whigh*stheta(kgrid)
          odifua(4,k,ista)=obsua(4,k,ista)-oanxua(4,k,ista)
          qbthe=wqlow*qback(4,ktab-1)+wqhigh*qback(4,ktab)
          theua(k,ista)=( (obsua(4,k,ista)/qobsua(4,k,ista))+           &
                       (oanxua(4,k,ista)/qbthe) ) /                     &
                       ((1./qobsua(4,k,ista)) + (1./qbthe))
          qobsua(4,k,ista)=qobsua(4,k,ista)/qbthe
        ELSE
          theua(k,ista)=wlow*stheta(kgrid-1)+whigh*stheta(kgrid)
        END IF
        IF(qualua(5,k,ista) > 0) THEN
          oanxua(5,k,ista)=wlow*sqv(kgrid-1)+whigh*sqv(kgrid)
          odifua(5,k,ista)=obsua(5,k,ista)-oanxua(5,k,ista)
          qobsua(5,k,ista)=qobsua(5,k,ista)/                            &
                 (wqlow*qback(5,ktab-1)+wqhigh*qback(5,ktab))
        END IF
! bad data for it is far below model ground level
        IF(whigh < -0.1 ) THEN
          DO ivar=1,nvar
            qualua(ivar,k,ista) = -999
          ENDDO
        ENDIF
!
      END DO

! shift ista to kntdom by mhu
!      IF( 1== 1 ) THEN
!      stnua(kntdom)=stnua(ista)
!      isrcua(kntdom)=isrcua(ista)
!      elevua(kntdom)=elevua(ista)
!      xua(kntdom)=xua(ista)
!      yua(kntdom)=yua(ista)
!      nlevsua(kntdom)=nlevsua(ista)
!      DO k=1,nlevsua(ista)
!        hgtua(k,kntdom)=hgtua(k,ista)
!        theua(k,kntdom)=theua(k,ista)
!        DO ivar=1,nvar
!          obsua(ivar,k,kntdom)=obsua(ivar,k,ista)
!          qobsua(ivar,k,kntdom)=qobsua(ivar,k,ista)
!          qualua(ivar,k,kntdom)=qualua(ivar,k,ista)
!          odifua(ivar,k,kntdom)=odifua(ivar,k,ista)
!          oanxua(ivar,k,kntdom)=oanxua(ivar,k,ista)
!        ENDDO
!      ENDDO
!     ENDIF

!
    ELSE
!
!      write(6,'(a,a)') stnua(ista),
!    :            ' outside ARPS domain, skipping'

!
!   The rule for "is an observation to be used" is if the ob is in the domain
!   somewhere.  For MPI, an ob that we might need might be outside of the local
!   processor domain, but inside of the large domain.  Make sure it stays
!   valid.
!

!      ireturn = -1
!      IF (mp_opt > 0 ) THEN
!          CALL findlc(nxlg,nylg,xslg,yslg,xua(ista),yua(ista),          &
!                ipt,jpt,ireturn)
!      END IF
! 
!      IF ( ireturn < 0 ) THEN
      IF ( indexua(ista) < 0 ) THEN  ! no one owns it
        isrcua(ista)=0

!     write(6,'(a,a,i4)') ' Extrapolation status returned from',
!    :                        ' subroutine findlc: ',ireturn
!      IF(ireturn.eq.-1) write(6,'(a)')
!    :      ' Extrapolation in x detected'
!      IF(ireturn.eq.-2) write(6,'(a)')
!    :      ' Extrapolation in y detected'
!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6/)')
!    :  ' location x= ',(0.001*xua(ista)),
!    :         '   y= ',(0.001*yua(ista)),
!    :  ' nearest pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))
!
         DO k=1,nzua
           DO ivar=1,nvar
             qualua(ivar,k,ista)=-9
           END DO
         END DO
       END IF

    END IF   ! inside domain?

  END DO

! wdt
  IF (mp_opt == 0) THEN
    WRITE(6,'(1x,i5,a,i5,a)') kntdom,' of ',nobsua,                     &
          ' upper air profiles are inside the ARPS domain'
  ELSE
    CALL flush(6)
    CALL mpbarrier
    DO k=0,nprocs-1
      IF (myproc == k .AND. kntdom > 0) THEN
        WRITE(6,'(1x,i5,a,i5,a,a,1x,i4)') kntdom,' of ',nobsua,         &
          ' upper air profiles are inside the ARPS domain',             &
          ' for proc ',myproc
        CALL flush(6)
      END IF
      CALL mpbarrier
    END DO
  END IF

!
! added by mhu
!
!  IF( 1==1 ) THEN
!   nobsua=kntdom
!  PRINT *, kntdom,' of ',nobsua,                    &
!           ' upper air profiles are inside the ARPS domain'
!  END IF

END SUBROUTINE grdtoua
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRDTORET                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grdtoret(nx,ny,nz,nztab,                                     &
           nzret,mxret,mxcolret,nvar,ncolret,                           &
           xs,ys,zp,anx,qback,hqback,nlvqback,                          &
           shght,su,sv,spres,stheta,sqv,                                &
           stnret,iret,xretc,yretc,hgtretc,                             &
           obsret,qobsret,qualret,nlevret,                              &
           odifret,oanxret,theretc,trnretc)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Gridded data are interpolated to retrieval column locations to
!  produce observation differences.  For now, bilinear interpolation
!  is used to create a vertical sounding at the obs location
!  and that is interpolated in the vertical to obtain the
!  bacground value.
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, July, 1995
!
!  MODIFICATION HISTORY:
!  Sept, 1996 (KB)
!  Added count-reporting diagnostic output.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz,nztab
  INTEGER :: nzret,mxret,mxcolret,nvar
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
!-----------------------------------------------------------------------
!
!  Background fields
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
  REAL :: qback(nvar,nztab)
  REAL :: hqback(nztab)
  INTEGER :: nlvqback
!
  REAL :: shght(nz)
  REAL :: su(nz)
  REAL :: sv(nz)
  REAL :: spres(nz)
  REAL :: stheta(nz)
  REAL :: sqv(nz)
!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnret(mxret)
  INTEGER :: iret(mxcolret)
  REAL :: xretc(mxcolret),yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: qobsret(nvar,nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: oanxret(nvar,nzret,mxcolret)
  REAL :: theretc(nzret,mxcolret)
  REAL :: trnretc(mxcolret)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,k,kgrid,ktab,kntdom
  INTEGER :: ipt,jpt,ireturn,nlevs
  REAL :: selev,hgtgnd,hgtagl
  REAL :: whigh,wlow,wqhigh,wqlow,qbthe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  kntdom=0
  DO ista=1,ncolret
!
!-----------------------------------------------------------------------
!
!  Find location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL findlc(nx,ny,xs,ys,xretc(ista),yretc(ista),                    &
                ipt,jpt,ireturn)
!
!-----------------------------------------------------------------------
!
!  Interpolate if extrapoltion is not indicated.
!
!-----------------------------------------------------------------------
!
    IF(ireturn >= 0) THEN
!      write(6,'(a,i4,a)')  'Column ',ista,
!    :            ' inside ARPS domain, processing'

      kntdom=kntdom+1

!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6)')
!    :  ' location x= ',(0.001*xua(ista)),
!    :         '   y= ',(0.001*yua(ista)),
!    :  ' found near pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2,/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))
!
!-----------------------------------------------------------------------
!
!  Interpolate (in the horizontal) for the whole vertical column.
!
!-----------------------------------------------------------------------
!
      CALL colint(nx,ny,nz,nvar,                                        &
             xs,ys,zp,xretc(ista),yretc(ista),ipt,jpt,anx,              &
             su,sv,stheta,spres,shght,sqv,selev,                        &
             nlevs)
!
!-----------------------------------------------------------------------
!
!  Convert pressure to ln(p) for interpolation
!
!-----------------------------------------------------------------------
!
      DO k=1,nlevs
        spres(k)=LOG(spres(k))
      END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate in the vertical at each obs level and compute
!  obs differences.
!
!-----------------------------------------------------------------------
!
      hgtgnd=0.5*(shght(2)+shght(1))
      trnretc(ista)=hgtgnd
      DO k=1,nlevret(ista)
        DO kgrid=2,nlevs-1
          IF(shght(kgrid) > hgtretc(k,ista)) GO TO 260
        END DO
        DO ivar=1,nvar
          qualret(ivar,k,ista)=-9
        END DO
        CYCLE
        260       CONTINUE
        whigh=(hgtretc(k,ista)-shght(kgrid-1))/                         &
               (shght(kgrid)-shght(kgrid-1))
        wlow=1.-whigh
!
        hgtagl=AMIN1(hqback(nlvqback),                                  &
                     AMAX1(hqback(1),(hgtretc(k,ista)-hgtgnd)))
        DO ktab=2,nlvqback-1
          IF(hqback(ktab) > hgtagl) EXIT
        END DO
!        281       CONTINUE
        wqhigh=(hgtagl-hqback(ktab-1))/                                 &
               (hqback(ktab)-hqback(ktab-1))
        wqlow=1.-wqhigh
!
        IF(qualret(1,k,ista) > 0) THEN
          oanxret(1,k,ista)=wlow*su(kgrid-1)+whigh*su(kgrid)
          odifret(1,k,ista)=obsret(1,k,ista)-oanxret(1,k,ista)
          qobsret(1,k,ista)=qobsret(1,k,ista)/                          &
                 (wqlow*qback(1,ktab-1)+wqhigh*qback(1,ktab))
        END IF
        IF(qualret(2,k,ista) > 0) THEN
          oanxret(2,k,ista)=wlow*sv(kgrid-1)+whigh*sv(kgrid)
          odifret(2,k,ista)=obsret(2,k,ista)-oanxret(2,k,ista)
          qobsret(2,k,ista)=qobsret(2,k,ista)/                          &
                 (wqlow*qback(2,ktab-1)+wqhigh*qback(2,ktab))
        END IF
        IF(qualret(3,k,ista) > 0) THEN
          oanxret(3,k,ista)=EXP( wlow*spres(kgrid-1)+                   &
                                 whigh*spres(kgrid) )
          odifret(3,k,ista)=obsret(3,k,ista)-oanxret(3,k,ista)
          qobsret(3,k,ista)=qobsret(3,k,ista)/                          &
                 (wqlow*qback(3,ktab-1)+wqhigh*qback(3,ktab))
        END IF
        IF(qualret(4,k,ista) > 0) THEN
          oanxret(4,k,ista)=wlow*stheta(kgrid-1)+whigh*stheta(kgrid)
          odifret(4,k,ista)=obsret(4,k,ista)-oanxret(4,k,ista)
          qbthe=wqlow*qback(4,ktab-1)+wqhigh*qback(4,ktab)
          theretc(k,ista)=( (obsret(4,k,ista)/qobsret(4,k,ista))+       &
                       (oanxret(4,k,ista)/qbthe) ) /                    &
                       ((1./qobsret(4,k,ista)) + (1./qbthe))
          qobsret(4,k,ista)=qobsret(4,k,ista)/qbthe
        ELSE
          theretc(k,ista)=wlow*stheta(kgrid-1)+whigh*stheta(kgrid)
        END IF
        IF(qualret(5,k,ista) > 0) THEN
          oanxret(5,k,ista)=wlow*sqv(kgrid-1)+whigh*sqv(kgrid)
          odifret(5,k,ista)=obsret(5,k,ista)-oanxret(5,k,ista)
          qobsret(5,k,ista)=qobsret(5,k,ista)/                          &
                 (wqlow*qback(5,ktab-1)+wqhigh*qback(5,ktab))
        END IF
      END DO
!
!       
! shift ista to kntdom by mhu
!      stnret(kntdom)=stnret(ista)
!      iret(kntdom)=iret(ista)
!      nlevret(kntdom)=nlevret(ista)
!      xretc(kntdom)=xretc(ista)
!      yretc(kntdom)=yretc(ista)
!      DO k=1,nlevret(kntdom)
!        hgtretc(k,kntdom)=hgtretc(k,ista)
!        theretc(k,kntdom)=theretc(k,ista)
!        DO ivar=1,nvar
!          obsret(ivar,k,kntdom)=obsret(ivar,k,ista)                     
!          qobsret(ivar,k,kntdom)=qobsret(ivar,k,ista)
!          qualret(ivar,k,kntdom)=qualret(ivar,k,ista)
!          odifret(ivar,k,kntdom)=odifret(ivar,k,ista)
!          oanxret(ivar,k,kntdom)=oanxret(ivar,k,ista)
!        ENDDO
!      ENDDO
!
    ELSE
!
!      write(6,'(a,a)') stnret(ista),
!    :            ' outside ARPS domain, skipping'

      iret(ista)=0

!     write(6,'(a,a,i4)') ' Extrapolation status returned from',
!    :                        ' subroutine findlc: ',ireturn
!      IF(ireturn.eq.-1) write(6,'(a)')
!    :      ' Extrapolation in x detected'
!      IF(ireturn.eq.-2) write(6,'(a)')
!    :      ' Extrapolation in y detected'
!      write (6,'(2x,a,f12.2,a,f12.2,/
!    :             2x,a,i6,a,i6/)')
!    :  ' location x= ',(0.001*xretc(ista)),
!    :         '   y= ',(0.001*yretc(ista)),
!    :  ' nearest pt i= ',ipt,'   j= ',jpt
!      write (6,'(12x,a,f12.2,a,f12.2/)')
!    :  'x= ',(0.001*xs(ipt)),'   y= ',(0.001*ys(jpt))
!
      DO k=1,nzret
        DO ivar=1,nvar
          qualret(ivar,k,ista)=-9
        END DO
      END DO

    END IF   ! inside domain?

  END DO
  IF (kntdom > 0) THEN
    WRITE(6,'(//,1x,i5,a,i5,a)') kntdom,' of ',ncolret,                 &
          ' retrieval columns are inside the ARPS domain'
  END IF
!  ncolret=kntdom
!  write(*,*) ' Use retrival data in the analysis domain',ncolret

  RETURN
END SUBROUTINE grdtoret
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE COLINT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE colint(nx,ny,nz,nvar,                                        &
           xs,ys,zp,xpt,ypt,ipt,jpt,anx,                                &
           su,sv,stheta,spres,shght,sqv,selev,                          &
           nlevs)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolates ARPS history data in the horizontal to create
!  a column of data located at point xpt, ypt.
!
!  Bilinear interpolation is used.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  October, 1992 (K. Brewster)
!  Conversion to ARPS 3.0.
!
!  October, 1994 (K. Brewster)
!  Conversion to ARPS 4.0.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx,ny,nz Dimensions of ARPS grids.
!
!    xs       x coordinate of scalar points in physical/comp. space (m)
!    ys       y coordinate of scalar points in physical/comp. space (m)
!    zp       z coordinate of scalar grid points in physical space (m)
!
!    xpt      x coordinate of desired sounding (m)
!    ypt      y coordinate of desired sounding (m)
!
!    ipt      i index of grid point just west of xpt,ypt
!    jpt      j index of grid point just south of xpt,ypt
!
!    anx      Background field
!
!  OUTPUT:
!
!    su       Interpolated u wind component.  (m/s)
!    sv       Interpolated v wind component.  (m/s)
!    stheta   Interpolated potential temperature (K).
!    spres    Interpolated pressure. (Pascals)
!    shght    Interpolated height (meters)
!    sqv      Interpolated specific humidity (kg/kg)
!    selev    Interpolated surface elevation (m)
!    nlevs    Number of above-ground sounding levels.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
!  Arguments -- location data
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nvar     ! Dimensions of ARPS grids.
  REAL :: xs(nx)               ! x coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: ys(ny)               ! y coordinate of grid points in
                               ! physical/comp. space (m)
  REAL :: zp(nx,ny,nz)         ! z coordinate of grid points in
                               ! physical space (m)
  REAL :: xpt                  ! location to find in x coordinate (m)
  REAL :: ypt                  ! location to find in y coordinate (m)
  INTEGER :: ipt               ! i index to the west of desired
                               ! location
  INTEGER :: jpt               ! j index to the south of desired
                               ! location
!
!-----------------------------------------------------------------------
!
!  Arguments -- background field
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
!
!-----------------------------------------------------------------------
!
!  Arguments -- Extracted sounding variables
!
!-----------------------------------------------------------------------
!
  REAL :: su(nz),sv(nz),stheta(nz),sqv(nz)
  REAL :: spres(nz),shght(nz)
  REAL :: selev
  INTEGER :: nlevs
!
!-----------------------------------------------------------------------
!
!  Functions called
!
!-----------------------------------------------------------------------
!
  REAL :: aint2d
!
!-----------------------------------------------------------------------
!
!  Include files
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
  INTEGER :: k,in,jn
  REAL :: delx,ddx,dely,ddy,w1,w2,w3,w4
  REAL :: t2,t3,hmid,tmid
!
!-----------------------------------------------------------------------
!
!  Find corner weights
!
!-----------------------------------------------------------------------
!
  in=ipt+1
  delx=xs(in)-xs(ipt)
  IF(ABS(delx) > 0.) THEN
    ddx=(xpt-xs(ipt))/delx
  ELSE
    ddx=0.
  END IF

  jn=jpt+1
  dely=ys(jn)-ys(jpt)
  IF(ABS(dely) > 0.) THEN
    ddy=(ypt-ys(jpt))/dely
  ELSE
    ddy=0.
  END IF

  w1=(1.-ddx)*(1.-ddy)
  w2=ddx*(1.-ddy)
  w3=ddx*ddy
  w4=(1.-ddx)*ddy
!
!-----------------------------------------------------------------------
!
!  Interpolate all variables at all levels.
!
!-----------------------------------------------------------------------
!
  nlevs=nz-1
  DO k=1,nz-1
    shght(k)=                                                           &
        aint2d(nx,ny,nz,    zp,ipt,jpt,k,in,jn,w1,w2,w3,w4)
    su(k)=                                                              &
        aint2d(nx,ny,nz, anx(1,1,1,1),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sv(k)=                                                              &
        aint2d(nx,ny,nz, anx(1,1,1,2),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    spres(k)=                                                           &
        aint2d(nx,ny,nz, anx(1,1,1,3),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    stheta(k)=                                                          &
        aint2d(nx,ny,nz, anx(1,1,1,4),ipt,jpt,k,in,jn,w1,w2,w3,w4)
    sqv(k)=                                                             &
        aint2d(nx,ny,nz, anx(1,1,1,5),ipt,jpt,k,in,jn,w1,w2,w3,w4)
  END DO
!
!-----------------------------------------------------------------------
!
!  Get height at scalar points, since zp was defined at w points.
!
!-----------------------------------------------------------------------
!
  selev=shght(2)
  DO k=1,nz-1
    shght(k)=0.5*(shght(k+1)+shght(k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Get a value at the surface, by linearly interpolating
!  between the 1st and second levels.
!
!-----------------------------------------------------------------------
!
  w2=(selev-shght(1))/(shght(2)-shght(1))
  w1=1.-w2
  su(1)=w1*    su(1) + w2*    su(2)
  sv(1)=w1*    sv(1) + w2*    sv(2)
  stheta(1)=w1*stheta(1) + w2*stheta(2)
  sqv(1)=w1*   sqv(1) + w2*   sqv(2)
  shght(1)=selev
!
!-----------------------------------------------------------------------
!
!  Integrate downward to get the pressure at level 1.
!
!-----------------------------------------------------------------------
!
  t3=stheta(3)*(spres(3)/100000.)**rddcp
  t2=stheta(2)*(spres(2)/100000.)**rddcp
  hmid=0.5*(shght(2)+shght(1))
  tmid=t3+(((shght(3)-hmid)/(shght(3)-shght(2)))*(t2-t3))
  spres(1)=spres(2)*EXP(g*(shght(2)-shght(1))/(rd*tmid))
  RETURN
END SUBROUTINE colint
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE QCDIFF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE qcdiff(nvar,nvarrad,nvarradin,mxsng,nsrcsng,                 &
           indexsng,indexua,indexrad,indexret,                          &
           usesng,useua,userad,                                         &
           nzua,mxua,nsrcua,                                            &
           nzrad,mxrad,mxcolrad,nsrcrad,                                &
           nzret,mxret,mxcolret,nsrcret,                                &
           nobsng,nobsua,ncolrad,ncolret,nam_var,                       &
           stnsng,isrcsng,hgtsng,obsng,oanxsng,odifsng,                 &
           qcthrsng,qualsng,                                            &
           stnua,isrcua,hgtua,obsua,oanxua,odifua,                      &
           qcthrua,qualua,nlevsua,                                      &
           stnrad,irad,isrcrad,hgtradc,obsrad,odifrad,                  &
           qcthrrad,qualrad,nlevrad,                                    &
           stnret,iret,isrcret,hgtretc,                                 &
           obsret,oanxret,odifret,                                      &
           qcthrret,qualret,nlevret,                                    &
           obsknt,rejknt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Observation differences are compared to expected errors
!  and rejected if they exceed a threshold.
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, July, 1995
!
!  MODIFICATION HISTORY:
!
!  Nov., 2002, Ming Hu
!  Deleted part converting QC citeria from Qv to RH for we don't need
!  it. We have got their QC citeria according to analysis type used.
! 
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  Sizing variables
!
  INTEGER :: nvar,nvarrad,nvarradin,mxsng,nsrcsng
  INTEGER :: nzua,mxua,nsrcua
  INTEGER :: nzrad,mxrad,mxcolrad,nsrcrad
  INTEGER :: nzret,mxret,mxcolret,nsrcret
!
  INTEGER :: nobsng,nobsua,ncolrad,ncolret
  CHARACTER (LEN=6) :: nam_var(nvar)

  INTEGER, INTENT(IN) :: indexsng(mxsng)
  INTEGER, INTENT(IN) :: indexua (mxua)
  INTEGER, INTENT(IN) :: indexrad(mxcolrad)
  INTEGER :: indexret(mxcolret)     ! place holder only, not used

  LOGICAL, INTENT(IN) :: usesng(mxsng)
  LOGICAL, INTENT(IN) :: useua (mxua)
  LOGICAL, INTENT(IN) :: userad(mxcolrad)

!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnsng(mxsng)
  INTEGER :: isrcsng(mxsng)
  REAL :: obsng(nvar,mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: qcthrsng(nvar,nsrcsng)
  INTEGER :: qualsng(nvar,mxsng)
  REAL :: oanxsng(nvar,mxsng)
  REAL :: odifsng(nvar,mxsng)
!
  CHARACTER (LEN=5) :: stnua(mxua)
  INTEGER :: isrcua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: qcthrua(nvar,nsrcua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
!
  CHARACTER (LEN=5) :: stnrad(mxrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  REAL :: hgtradc(nzrad,mxcolrad)
  REAL :: obsrad(nvarradin,nzrad,mxcolrad)
  REAL :: odifrad(nvarrad,nzrad,mxcolrad)
  REAL :: qcthrrad(nvarradin,nsrcrad)
  INTEGER :: qualrad(nvarrad,nzrad,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
!
  CHARACTER (LEN=5) :: stnret(mxret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: oanxret(nvar,nzret,mxcolret)
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: qcthrret(nvar,nsrcret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)

  REAL :: obsknt(nvar)
  REAL :: rejknt(nvar)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,klev,jsrc
  REAL    :: pctrej,vrdiff
  REAL    :: tk,qvsat,qcthrqv
  INTEGER :: mvar
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mvar = nvar - 1
!
!-----------------------------------------------------------------------
!
!  Check single-level observation differences
!
!-----------------------------------------------------------------------
!
  DO ivar=1,nvar
    obsknt(ivar)=0.
    rejknt(ivar)=0.
  END DO
  IF (mp_opt > 0) THEN
    CALL flush(6)
    CALL mpbarrier
  END IF
  DO ista=1,nobsng
!    IF ( mp_opt == 0 .OR. (mp_opt > 0 .AND. indexsng(ista) == myproc )) THEN
    IF ( usesng(ista) ) THEN
!write(100+myproc,*) 'station: ',ista, isrcsng(ista),qualsng(4,ista),odifsng(4,ista),qcthrsng(4,isrcsng(ista))

        IF(isrcsng(ista) > 0) THEN
          DO ivar=1,mvar
            IF(qualsng(ivar,ista) > 0) THEN
              IF ( indexsng(ista) == myproc ) obsknt(ivar)=obsknt(ivar)+1.
!write(100+myproc,*) 'count obs: ',ista, ivar, obsknt(ivar)
              IF(ABS(odifsng(ivar,ista)) > qcthrsng(ivar,isrcsng(ista)) ) THEN
                IF (indexsng(ista) == myproc) THEN
                  WRITE(6,'(a,a,a,a)')                                  &
                      ' Discarding ',nam_var(ivar),                     &
                      ' at ',stnsng(ista)
                  WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')           &
                      ' Hgt= ',hgtsng(ista),' Ob= ',obsng(ivar,ista),   &
                      ' Diff= ',odifsng(ivar,ista),                     &
                      ' Thresh=',qcthrsng(ivar,isrcsng(ista))
                  CALL flush(6)
                END IF
!write(100+myproc,*) 'reject obs: ',ista
                qualsng(ivar,ista)=-199
                IF ( indexsng(ista) == myproc ) rejknt(ivar)=rejknt(ivar)+1.
              END IF
            END IF
          END DO

!
!  Don't know why this is commented out.  ADAS wants it to be there.  Maybe
!  3dvar doesn't need it.  I'll "ifdef" it to be sure...kwt.
!

         IF(qualsng(5,ista) > 0) THEN
           IF ( indexsng(ista) == myproc ) obsknt(5)=obsknt(5)+1.
           tk=oanxsng(4,ista)*((oanxsng(3,ista)/p0)**rddcp)
           qvsat=f_qvsat(oanxsng(3,ista),tk)
           qcthrqv=qcthrsng(5,isrcsng(ista))*qvsat
           IF(ABS(odifsng(5,ista)) > qcthrqv) THEN
!             IF (mp_opt == 0 .OR. indexsng(ista) == myproc) THEN
             IF (indexsng(ista) == myproc) THEN
               WRITE(6,'(a,a,a,a)')                                     &
                   ' Discarding ',nam_var(5),                           &
                   ' at ',stnsng(ista)
               WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')              &
                   ' Hgt= ',hgtsng(ista),' Ob= ',obsng(5,ista),         &
                   ' Diff= ',odifsng(5,ista),                           &
                   ' Thresh=',qcthrqv
               IF (mp_opt > 0) CALL flush(6)
             END IF
             qualsng(5,ista)=-199
             IF ( indexsng(ista) == myproc ) rejknt(5)=rejknt(5)+1.
           END IF
         END IF

    !
    !-----------------------------------------------------------------------
    !
    !  Synchronize wind component flags
    !
    !-----------------------------------------------------------------------
    !
          IF(qualsng(1,ista) == -199) qualsng(2,ista)=-199
          IF(qualsng(2,ista) == -199) qualsng(1,ista)=-199
    
        END IF
    END IF
    CALL mpbarrier
  END DO
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    DO ivar=1,nvar
      CALL mpsumr(obsknt(ivar),1)
      CALL mpsumr(rejknt(ivar),1)
    END DO
  END IF
  IF (myproc == 0) THEN
    WRITE(6,'(//a/a)') ' QC Statistics for single-level data',          &
                     ' Variable     Obs      Rejected    Percent'
    DO ivar=1,nvar
      pctrej=0.
      IF(obsknt(ivar) > 0.) pctrej=100.*rejknt(ivar)/obsknt(ivar)
        WRITE(6,'(2x,a,i10,i10,f8.1)')                                  &
          nam_var(ivar),nint(obsknt(ivar)),nint(rejknt(ivar)),pctrej
    END DO
    WRITE(6,'(/a)') ' '
  END IF
!
!-----------------------------------------------------------------------
!
!    Check upper-air observation differences
!
!-----------------------------------------------------------------------
!
  DO ivar=1,nvar
    obsknt(ivar)=0.
    rejknt(ivar)=0.
  END DO
  IF (mp_opt > 0) THEN
    CALL flush(6)
    CALL mpbarrier
  END IF
  DO ista=1,nobsua
!    IF ( mp_opt == 0 .OR. (mp_opt > 0 .AND. indexua(ista) == myproc )) THEN
    IF ( useua(ista) ) THEN
        IF(isrcua(ista) > 0) THEN
          DO klev=1,nlevsua(ista)
            DO ivar=1,mvar
              IF(qualua(ivar,klev,ista) > 0) THEN
                IF (indexua(ista) == myproc) obsknt(ivar)=obsknt(ivar)+1.
                IF (ABS(odifua(ivar,klev,ista))                         &
                    > qcthrua(ivar,isrcua(ista))) THEN
                  IF (indexua(ista) == myproc) THEN
                    WRITE(6,'(a,a,a,a,a,i4)') ' Discarding UA ',        &
                      nam_var(ivar),' at sta ',stnua(ista),             &
                      ' klev =',klev
                    WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')         &
                      ' Hgt= ',hgtua(klev,ista),                        &
                      ' Ob= ',obsua(ivar,klev,ista),                    &
                      ' Diff= ',odifua(ivar,klev,ista),                 &
                      ' Thresh=',qcthrua(ivar,isrcua(ista))
                    CALL flush(6)
                  END IF
                  qualua(ivar,klev,ista)=-199
                  IF (indexua(ista) == myproc) rejknt(ivar)=rejknt(ivar)+1.
                END IF
              END IF
            END DO
!           IF(qualua(5,klev,ista) > 0) THEN
!             obsknt(5)=obsknt(5)+1.
!             tk=oanxua(4,klev,ista)*                                       &
!                      ((oanxua(3,klev,ista)/p0)**rddcp)
!             qvsat=f_qvsat(oanxua(3,klev,ista),tk)
!             qcthrqv=qcthrua(5,isrcua(ista))*qvsat
!             IF(ABS(odifua(5,klev,ista)) > qcthrqv) THEN
!               IF (mp_opt == 0 .OR. indexua(ista) == myproc)  THEN
!                 WRITE(6,'(a,a,a,a)')                                      &
!                   ' Discarding ',nam_var(5),                              &
!                   ' at ',stnua(ista)
!                 WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')               &
!                   ' Hgt= ',hgtua(klev,ista),                              &
!                   ' Ob= ',obsua(5,klev,ista),                             &
!                   ' Diff= ',odifua(5,klev,ista),                          &
!                   ' Thresh=',qcthrqv
!                 IF (mp_opt > 0) CALL flush(6)
!               END IF
!               qualua(5,klev,ista)=-199
!               rejknt(5)=rejknt(5)+1.
!             END IF
!           END IF
!
!-----------------------------------------------------------------------
!
!  Synchronize wind component flags
!
!-----------------------------------------------------------------------
!
            IF(qualua(1,klev,ista) == -199) qualua(2,klev,ista)=-199
            IF(qualua(2,klev,ista) == -199) qualua(1,klev,ista)=-199
    
          END DO
        END IF
    END IF
    CALL mpbarrier
  END DO
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    DO ivar=1,nvar
      CALL mpsumr(obsknt(ivar),1)
      CALL mpsumr(rejknt(ivar),1)
    END DO
  END IF
  IF (myproc == 0) THEN
    WRITE(6,'(//a/a)') ' QC Statistics for UA data',                    &
                     ' Variable     Obs      Rejected    Percent'
    DO ivar=1,nvar
      pctrej=0.
      IF(obsknt(ivar) > 0.) pctrej=100.*rejknt(ivar)/obsknt(ivar)
      WRITE(6,'(2x,a,i10,i10,f8.1)')                                    &
        nam_var(ivar),nint(obsknt(ivar)),nint(rejknt(ivar)),pctrej
    END DO
    WRITE(6,'(/a)') ' '
  END IF
!
!-----------------------------------------------------------------------
!
!    Check radar data observation differences
!    Here radial wind difference is computed and compared only.
!
!-----------------------------------------------------------------------
!
  obsknt(1)=0.
  rejknt(1)=0.

  CALL flush(6)
  CALL mpbarrier

  DO ista=1,ncolrad
!    IF ( mp_opt == 0 .OR. (mp_opt > 0 .AND. indexrad(ista) == myproc )) THEN
    IF ( userad(ista) ) THEN
      IF(irad(ista) > 0) THEN
        jsrc=isrcrad(irad(ista))
        DO klev=1,nlevrad(ista)
        IF(qualrad(1,klev,ista) > 0 .AND. qualrad(2,klev,ista) > 0 ) THEN
          IF ( indexrad(ista) == myproc ) obsknt(1)=obsknt(1)+1.
          vrdiff=SQRT(                                                  &
                 odifrad(1,klev,ista)*odifrad(1,klev,ista)+             &
                 odifrad(2,klev,ista)*odifrad(2,klev,ista))

          IF( ABS(vrdiff) > qcthrrad(2,jsrc) ) THEN
            IF ( indexrad(ista) == myproc ) THEN
              WRITE(6,'(a,a,i5,a,i4)') ' Discarding radar Vr',            &
                ' at col ',ista,' klev =',klev
              WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')                 &
                ' Hgt= ',hgtradc(klev,ista),                              &
                ' Ob= ',obsrad(2,klev,ista),                              &
                ' Diff= ',vrdiff,                                         &
                ' Thresh=',qcthrrad(2,jsrc)
              CALL flush(6)
            END IF
            qualrad(1,klev,ista)=-199
            qualrad(2,klev,ista)=-199
            IF ( indexrad(ista) == myproc ) rejknt(1)=rejknt(1)+1.
          END IF
        END IF
        END DO
      END IF
    END IF
!
!   Can't "mpbarrier" here, as each radar only has the obs that it needs.
!   Some radars will have ncolrad=0, so won't get here.  Since I've not
!   seen velocities thresholded out, I don't think I have to worry about
!   many procs writing messages at the same time.
!
!   CALL mpbarrier
  END DO

  IF (mp_opt > 0) THEN
    DO ivar=1,nvar
      CALL mpsumr(obsknt(ivar),1)
      CALL mpsumr(rejknt(ivar),1)
    END DO
  END IF

  IF (myproc == 0) THEN
    WRITE(6,'(//a/a)') ' QC Statistics for radar data',                 &
                     ' Variable     Obs      Rejected    Percent'
    pctrej=0.
    IF(obsknt(1) > 0.) pctrej=100.*rejknt(1)/obsknt(1)
    WRITE(6,'(2x,a,4x,i10,i10,f8.1)')                                    &
          'Vr',nint(obsknt(1)),nint(rejknt(1)),pctrej
    WRITE(6,'(/a)') ' '
  END IF
!
!-----------------------------------------------------------------------
!
!    Check retrieval observation differences
!
!-----------------------------------------------------------------------
!
  DO ivar=1,nvar
    obsknt(ivar)=0.
    rejknt(ivar)=0.
  END DO
  IF (mp_opt > 0) THEN
    CALL flush(6)
    CALL mpbarrier
  END IF
  DO ista=1,ncolret
    IF(iret(ista) > 0) THEN
      jsrc=isrcret(iret(ista))
      DO klev=1,nlevret(ista)
        DO ivar=1,nvar-1
          IF(qualret(ivar,klev,ista) > 0) THEN
            obsknt(ivar)=obsknt(ivar)+1.
            IF(ABS(odifret(ivar,klev,ista)) >   qcthrret(ivar,jsrc) ) THEN
!
! RED FLAG: no "indexret" variable, so no order!!!
!
              WRITE(6,'(a,a,a,i5,a,i4)') ' Discarding retrieval ',      &
                  nam_var(ivar),' at col ',ista,' klev =',klev
              WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')               &
                  ' Hgt= ',hgtretc(klev,ista),                          &
                  ' Ob= ',obsret(ivar,klev,ista),                       &
                  ' Diff= ',odifret(ivar,klev,ista),                    &
                  ' Thresh=',qcthrret(ivar,jsrc)
              qualret(ivar,klev,ista)=-199
              rejknt(ivar)=rejknt(ivar)+1.
            END IF
          END IF
        END DO
        IF(qualret(5,klev,ista) > 0) THEN
          obsknt(5)=obsknt(5)+1.
!mhu      tk=oanxret(4,klev,ista)*                                      &
!mhu                ((oanxret(3,klev,ista)/p0)**rddcp)
!mhu      qvsat=f_qvsat(oanxret(3,klev,ista),tk)
!mhu      qcthrqv=qcthrret(5,jsrc)*qvsat
          IF(ABS(odifret(5,klev,ista)) > qcthrqv) THEN
            IF (myproc == 0)                                            &
            WRITE(6,'(a,a,a,i5,a,i4)') ' Discarding retrieval ',        &
                  nam_var(5),' at col ',ista,' klev =',klev
            IF (myproc == 0)                                            &
            WRITE(6,'(a,f8.0,a,f16.6,a,f16.6,a,f16.6)')                 &
                ' Hgt= ',hgtretc(klev,ista),                            &
                ' Ob= ',obsret(5,klev,ista),                            &
                ' Diff= ',odifret(5,klev,ista),                         &
                ' Thresh=',qcthrqv
            qualret(5,klev,ista)=-199
            rejknt(5)=rejknt(5)+1.
          END IF
        END IF
!
!-----------------------------------------------------------------------
!
!  Synchronize wind component flags
!
!-----------------------------------------------------------------------
!
        IF(qualret(1,klev,ista) == -199) qualret(2,klev,ista)=-199
        IF(qualret(2,klev,ista) == -199) qualret(1,klev,ista)=-199

      END DO
    END IF
    CALL mpbarrier
  END DO
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    DO ivar=1,nvar
      CALL mpsumr(obsknt(ivar),1)
      CALL mpsumr(rejknt(ivar),1)
    END DO
  END IF
  IF (myproc == 0) THEN
    WRITE(6,'(//a/a)') ' QC Statistics for retrieval data',             &
                     ' Variable     Obs      Rejected    Percent'
    DO ivar=1,nvar
      pctrej=0.
      IF(obsknt(ivar) > 0.) pctrej=100.*rejknt(ivar)/obsknt(ivar)
      WRITE(6,'(2x,a,i10,i10,f8.1)')                                    &
        nam_var(ivar),nint(obsknt(ivar)),nint(rejknt(ivar)),pctrej
    END DO
    WRITE(6,'(/a)') ' '
  END IF
  RETURN
END SUBROUTINE qcdiff
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRDTOSNG                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE qv2rhobs(nvar,mxsng,nobsng,obsng,qualsng,                    &
                    nzua,mxua,nobsua,obsua,qualua,nlevsua,              &
                    nzret,mxcolret,ncolret,obsret,qualret,nlevret)      
!                   
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  
!   convert qv to rh  in observation array
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Ming Hu, CAPS, Nov, 2002
!
!  MODIFICATION HISTORY:
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: mxsng,nvar,nobsng
  INTEGER :: nzua,mxua,nobsua
  INTEGER :: nzret,mxcolret,ncolret
!
!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  REAL :: obsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)

  REAL :: obsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)

  REAL :: obsret(nvar,nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
!
  include 'phycst.inc'
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ivar,k
  REAL :: f_qvsat, temp_qvsat,temp
!
     print*,'qobsng=',nobsng
!     write(*,'(5f15.7)') (qobsng(4,ista),ista=1,nobsng)
!-----------------------------------------------------------------------
!
!  Beginning of executable code
!
!-----------------------------------------------------------------------
!
!  for single level data
!
  DO ista=1,nobsng
    IF(qualsng(5,ista) > 0 .and.                                 &
       qualsng(4,ista) > 0 .and.                                 &
       qualsng(3,ista) > 0 ) THEN
      temp=obsng(4,ista)*(obsng(3,ista)/p0)**rddcp
      temp_qvsat=f_qvsat(obsng(3,ista),temp)
      obsng(5,ista)=obsng(5,ista)/temp_qvsat
      IF( obsng(5,ista) > 1.0 ) obsng(5,ista)=1.0
      IF( obsng(5,ista) < 0.0 ) obsng(5,ista)=0.0
    ELSE
      qualsng(5,ista) = -999
    ENDIF
  END DO
!
! for upper level moisture observation
!
  DO ista=1,nobsua
    DO k=1,nlevsua(ista)
      IF(qualua(3,k,ista) > 0 .and.                             &
         qualua(4,k,ista) > 0 .and.                             &
         qualua(5,k,ista) > 0) THEN
        temp=obsua(4,k,ista)*(obsua(3,k,ista)/p0)**rddcp
        temp_qvsat=f_qvsat(obsua(3,k,ista),temp)
        obsua(5,k,ista)=obsua(5,k,ista)/temp_qvsat
        IF( obsua(5,k,ista) > 1.0 ) obsua(5,k,ista)=1.0
        IF( obsua(5,k,ista) < 0.0 ) obsua(5,k,ista)=0.0
      ELSE
        qualua(5,k,ista) = -999
      END IF
    ENDDO
  END DO
!
!  For retrieval moisture
!
  DO ista=1,ncolret
    DO k=1,nlevret(ista)
      IF(qualret(3,k,ista) > 0 .and.                          &
         qualret(4,k,ista) > 0 .and.                          &
         qualret(5,k,ista) > 0) THEN
        temp=obsret(4,k,ista)*(obsret(3,k,ista)/p0)**rddcp
        temp_qvsat=f_qvsat(obsret(3,k,ista),temp)
        obsret(5,k,ista)=obsret(5,k,ista)/temp_qvsat
        IF( obsret(5,k,ista) > 1.0 ) obsret(5,k,ista)=1.0
        IF( obsret(5,k,ista) < 0.0 ) obsret(5,k,ista)=0.0
      ELSE
        qualret(5,k,ista)=-999
      ENDIF
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE qv2rhobs


SUBROUTINE createua(nx,ny,nz,nzua,mxua,nvar,nobsua,xs,ys,zp,       &
                    xua,yua,hgtua,qobsua,qualua,nlevsua,odifua)
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRDTOUA                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!  PURPOSE:
!  create my own up level data
!
!  AUTHOR:
!
!  Ming Hu, CAPS, March, 2003
! 
!  MODIFICATION HISTORY:
! 
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: nzua,mxua,nvar
  INTEGER :: nobsua
!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
!
!-----------------------------------------------------------------------
!
!  Station variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnua(mxua)
  INTEGER :: isrcua(mxua)
  REAL :: elevua(mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  REAL :: odifua(nvar,nzua,mxua)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  REAL, DIMENSION(:,:,:), allocatable ::  uvw1
  REAL, DIMENSION(:,:,:), allocatable ::  uvw2

  INTEGER :: ista,ivar,k,kntdom
  INTEGER :: i,j,ibegin,jbegin,is,js
  INTEGER :: nwx,nwy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nwx=100
  nwy=100
  allocate( uvw1(nwx,nwy,2) )
  allocate( uvw2(nwx,nwy,2) )
    uvw1=0
    uvw2=0
    OPEN(13,file='../wind/psiwind.dat',form='unformatted')
     READ(13) uvw1
     READ(13) uvw2
     READ(13)
    CLOSE(13)

   ibegin=40
   jbegin=20
   kntdom=0
   ista=0
   DO j=1,nwy
   DO i=1,nwx
     ista=ista+1
     if( ista > mxua ) THEN
       write(*,*) 'too many ua'
       stop  234
     ENDIF
     is=ibegin+i
     js=jbegin+j
     nlevsua(ista)=1
     xua(ista)=xs(is)
     yua(ista)=ys(js)
     DO k=1,nlevsua(ista)
        hgtua(k,ista)=(zp(is,js,nz/2)+zp(is,js,nz/2+1))/2.0
        DO ivar=1,nvar
          qualua(ivar,k,ista)=-999
        enddo
        qobsua(1,k,ista)=0.5
        qualua(1,k,ista)=100
        odifua(1,k,ista)=uvw1(i,j,1)
!        qobsua(2,k,ista)=0.5
!        qualua(2,k,ista)=100
!        odifua(2,k,ista)=uvw2(i,j,1)
     ENDDO
   ENDDO
   ENDDO

   nobsua=ista
   IF(1 == 1 ) THEN
        xua(1)=(xs(107)+xs(108))/2.0
        yua(1)=ys(87)
        xua(3)=(xs(109)+xs(110))/2.0
        yua(3)=ys(87)
        xua(2)=(xs(107)+xs(108))/2.0
        yua(2)=ys(89)
        xua(4)=(xs(109)+xs(110))/2.0
        yua(4)=ys(89)
        hgtua(1,1)=(zp(108,87,26)+zp(108,87,27))/2.0
        hgtua(1,2)=(zp(108,89,26)+zp(108,89,27))/2.0
        hgtua(1,3)=(zp(109,87,26)+zp(109,87,27))/2.0
        hgtua(1,4)=(zp(110,89,26)+zp(110,89,27))/2.0
        qobsua(1,1,1)=0.5
        qualua(1,1,1)=100
        odifua(1,1,1)=10.0
        qobsua(1,1,2)=0.5
        qualua(1,1,2)=100
        odifua(1,1,2)=-10.0
        qobsua(1,1,3)=0.5
        qualua(1,1,3)=100
        odifua(1,1,3)=-10.0
        qobsua(1,1,4)=0.5
        qualua(1,1,4)=100
        odifua(1,1,4)=10.0
        nobsua=4
   ENDIF


END
