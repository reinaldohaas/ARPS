!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PREPRETR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE prepretr(nx,ny,nz,nvar,                                      &
           nzret,mxret,mxcolret,mxztab,nsrcret,                         &
           nretfil,fretname,                                            &
           isrcret,srcret,nlvrttab,qsrcret,hrtqsrc,                     &
           stnret,latret,lonret,elvret,                                 &
           latretc,lonretc,xretc,yretc,iret,nlevret,                    &
           hgtretc,obsret,qrret,qobsret,qualret,                        &
           rmiss,ncolret,tem1,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Control reading columns of data from the ARPS temperaturea
!  and velocity retrieval routines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  April, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nvar
  INTEGER :: nzret,mxret,mxcolret,mxztab,nsrcret

  INTEGER :: nretfil
  CHARACTER (LEN=256) :: fretname(nretfil)
  INTEGER :: isrcret(0:mxret)
  CHARACTER (LEN=8) :: srcret(mxret)
  INTEGER :: nlvrttab(nsrcret)
  REAL :: qsrcret(nvar,mxztab,nsrcret)
  REAL :: hrtqsrc(mxztab,nsrcret)
  CHARACTER (LEN=5) :: stnret(mxret)
  REAL :: latret(mxret)
  REAL :: lonret(mxret)
  REAL :: elvret(mxret)
  REAL :: latretc(mxcolret)
  REAL :: lonretc(mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: nlevret(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: qrret(nzret,mxcolret)
  REAL :: qobsret(nvar,nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  REAL :: rmiss
  INTEGER :: ncolret
  REAL :: tem1(nx,ny,nz)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ista,ilev,ivar,isrc,ktab
  REAL :: wthi,wtlo
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO ista=1,mxcolret
    iret(ista)=0
    DO ilev=1,nzret
      DO ivar=1,nvar
        obsret(ivar,ilev,ista)=rmiss
        qobsret(ivar,ilev,ista)=999999.
        qualret(ivar,ilev,ista)=0
      END DO
    END DO
  END DO
!
  CALL rdretcol(nx,ny,nz,nvar,                                          &
            mxret,nzret,mxcolret,nretfil,fretname,                      &
            srcret,isrcret(1),stnret,latret,lonret,elvret,              &
            latretc,lonretc,iret,nlevret,hgtretc,obsret,qrret,          &
            ncolret,istatus,tem1)
!
!-----------------------------------------------------------------------
!
!  Get x and y locations of each retrieval column
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(ncolret,1,latretc,lonretc,xretc,yretc)
!
!-----------------------------------------------------------------------
!
!  Set qobs based on source and height
!
!-----------------------------------------------------------------------
!
  DO ista=1,ncolret
    IF(iret(ista) > 0) THEN
      isrc=isrcret(iret(ista))
      DO ilev=1,nlevret(ista)
        DO ktab=2,nlvrttab(isrc)-1
          IF(hrtqsrc(ktab,isrc) > hgtretc(ilev,ista)) EXIT
        END DO
!        226       CONTINUE
        wthi=  (hgtretc(ilev,ista)-hrtqsrc(ktab-1,isrc))/               &
             (hrtqsrc(ktab,isrc)-hrtqsrc(ktab-1,isrc))
        wthi=AMAX1(wthi,0.0)
        wthi=AMIN1(wthi,1.0)
        wtlo=1.0-wthi
        DO ivar=1,nvar
          qobsret(ivar,ilev,ista)=                                      &
                 wthi*qsrcret(ivar,ktab,  isrc) +                       &
                 wtlo*qsrcret(ivar,ktab-1,isrc)
          IF(obsret(ivar,ilev,ista) > -999.) qualret(ivar,ilev,ista)=10
        END DO
      END DO
    END IF
  END DO
  RETURN
END SUBROUTINE prepretr
!

SUBROUTINE retmcro(nx,ny,nz,                                            &
           mx_ret,nsrc_ret,nvar_anx,nz_ret,mx_colret,                   &
           srcret,isrcret,iret,nlevret,                                 &
           xretc,yretc,hgtretc,qrret,ncolret,                           &
           dzfill,                                                      &
           xs,ys,zs,qr)
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: mx_ret,nsrc_ret,nvar_anx,nz_ret,mx_colret
!
!-----------------------------------------------------------------------
!
!  Retrieval observation variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcret(nsrc_ret)
  INTEGER :: isrcret(0:mx_ret)
  INTEGER :: iret(mx_colret)
  INTEGER :: nlevret(mx_colret)

  REAL :: xretc(mx_colret)
  REAL :: yretc(mx_colret)
  REAL :: hgtretc(nz_ret,mx_colret)
  REAL :: qrret(nz_ret,mx_colret)
!
  INTEGER :: ncolret
  REAL :: dzfill
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  REAL :: qr(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,icol,jret
  INTEGER :: klow,knex,klast
  REAL :: dx,dy,dstlim2,dist,qrlow,qrintr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dx=xs(2)-xs(1)
  dy=ys(2)-ys(1)
  dstlim2=0.95*(dx*dx+dy*dy)
!
!-----------------------------------------------------------------------
!
!  For each horizontal grid point, use retrieval columns within
!  the threshold distance defined by dstlim2 to set clouds and rain.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      DO icol=1,ncolret
        jret=iret(icol)
        IF(isrcret(jret) > 0) THEN
          dist=((xretc(icol)-xs(i))*(xretc(icol)-xs(i)) +               &
                (yretc(icol)-ys(j))*(yretc(icol)-ys(j)))
          IF(dist < dstlim2) THEN
            DO klow=1,nlevret(icol)
              IF(qrret(klow,icol) > 0.) GO TO 451
            END DO
            CYCLE
            451         CONTINUE
!
!-----------------------------------------------------------------------
!
!  Found a postive qr in retrieval data column.
!
!-----------------------------------------------------------------------
!
            480         CONTINUE
            DO k=1,nz-1
              IF(zs(i,j,k) >= hgtretc(klow,icol)) GO TO 501
            END DO
            CYCLE
            501         CONTINUE
            qr(i,j,k)=MAX(qr(i,j,k),qrret(klow,icol))
            qrlow=qrret(klow,icol)
            klast=k
!
!-----------------------------------------------------------------------
!
!  Find next positive qr in retrieval data column.
!
!-----------------------------------------------------------------------
!
            520         CONTINUE
            DO knex=klow+1,nlevret(icol)
              IF(qrret(klow,icol) > 0.) GO TO 551
            END DO
            CYCLE
            551         CONTINUE
!
!-----------------------------------------------------------------------
!
!  Found another positive qr, so determine which
!  ARPS grid points are affected.   Either
!  1) Fill from the previous reflectivity level  ..or..
!  2) Set a new klow and klast and treat as if this were
!     the first point for this column.
!
!-----------------------------------------------------------------------
!
            IF((hgtretc(knex,icol)-hgtretc(klow,icol)) <                &
                                                dzfill) THEN
              DO k=klast+1,nz-1
                IF(zs(i,j,k) <= hgtretc(knex,icol)) THEN
                  qrintr=(qrret(knex,icol)-qrlow)                       &
                      /(hgtretc(knex,icol)-hgtretc(klow,icol))          &
                      *(zs(i,j,k)-hgtretc(klow,icol))                   &
                      +qrlow
                  qr(i,j,k)=MAX(qr(i,j,k),qrintr)
                  klast=k
                ELSE
                  GO TO 601
                END IF
              END DO
              CYCLE
              601           CONTINUE
              qrlow=qrret(knex,icol)
              klow=knex
              GO TO 520
            ELSE
              qrlow=qrret(knex,icol)
              klow=knex
              GO TO 480
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE retmcro
