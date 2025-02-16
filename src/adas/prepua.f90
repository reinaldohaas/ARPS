!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PREPUAOBS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE prepuaobs(nx,ny,nz,nvar,nzua,mxua,mxztab,nsrcua,mxuafile,    &
           anx,xs,ys,zp,su,sv,st,spres,shght,sqv,sqvsat,                &
           nuafile,uafname,srcua,                                       &
           stnua,elevua,xua,yua,hgtua,obsua,                            &
           qsrcua,hgtqsrc,nlevtab,                                      &
           qobsua,qualua,isrcua,nlevsua,                                &
           rmiss,nobsua,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Control reading upper-air observations and preparation
!  of upper-air data for analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (K. Brewster)
!  Added full documentation.
!
!  Nov, 1997 (KB)
!  Changed arguments to be arrays of upper-air data filenames
!  rather than one each of a raob and profiler file.
!
!  July, 2011 (KB)
!  Moved the setting of the source into rdprof subroutines
!  This allows multiple profiler sources within the .pro files.
!
!  August, 2011 (KB)
!  Moved the setting of the source into rdraob subroutine.
!  This allows multiple sounding sources within the .snd files.
!
!-----------------------------------------------------------------------
!
!  nvar      Number of analysis variables
!  nzua      Maximum number of levels
!  mxua      Maximum number of multiple-level profiles
!  mxztab    Maximum number of levels in the data error table
!  nsrcua    Number of multiple-level data sources
!  mxuafile  Maximum number of multiple-level data files
!  nuafile   Number of multiple-level data files
!  uafname   File names of data files
!  stnua     Station name
!  elevua    Station elevation
!  xua       X-coordinate of station
!  yua       Y-coordinate of station
!  hgtua     Height of each level in data
!  obsua     Multiple-level observations
!  qsrcua    Standard error for each data source
!  hgtqsrc   Height of data in error table
!  nlevtab   Number of levels in error table
!  qobsua    Standard error for each datum
!  qualua    Quality indicator for each datum
!  isrcua    Source number for each station
!  nlevsua   Number of level for each station
!  rmiss     Missing data flag value
!  nobsua    Number of observations
!  istatus   Return status
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nzua,mxua,mxztab,nsrcua
  INTEGER :: mxuafile,nuafile

  REAL :: anx(nvar,nx,ny,nz)
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: su(nz)
  REAL :: sv(nz)
  REAL :: st(nz)
  REAL :: spres(nz)
  REAL :: shght(nz)
  REAL :: sqv(nz)
  REAL :: sqvsat(nz)

  CHARACTER (LEN=256) :: uafname(mxuafile)
  CHARACTER (LEN=8)   :: srcua(nsrcua)
  CHARACTER (LEN=5)   :: stnua(mxua)
  REAL :: elevua(mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: qsrcua(nvar,mxztab,nsrcua)
  REAL :: hgtqsrc(mxztab,nsrcua)
  REAL :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: nlevtab(nsrcua)
  INTEGER :: nobsua
  REAL :: rmiss
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=12) :: suffix
  CHARACTER (LEN=8) :: dfltsrc
  INTEGER :: ista,ilev,ivar,ifile,isrc,nprev,ntotal,ktab
  INTEGER :: maxsuf,lenfnm,dotloc,lensuf
  REAL :: wthi,wtlo
  LOGICAL :: found
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ntotal = 0

  DO ista=1,mxua
    isrcua(ista)=0
    DO ilev=1,nzua
      DO ivar=1,nvar
        qobsua(ivar,ilev,ista)=999999.
      END DO
    END DO
  END DO
!
  nprev=0
  maxsuf=LEN(suffix)
!
  DO ifile=1,nuafile
    lenfnm=LEN(uafname(ifile))
    CALL strlnth(uafname(ifile),lenfnm)
    CALL exsufx(uafname(ifile),lenfnm,suffix,maxsuf,dotloc,lensuf)

    IF(lensuf == 3 .AND. suffix(1:3) == 'snd') THEN

      dfltsrc = 'NWS RAOB'
      CALL rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 3 .AND. suffix(1:3) == 'pro') THEN

      dfltsrc = 'WPDN PRO'
      CALL rdprof(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 3 .AND. suffix(1:3) == 'vad') THEN

      dfltsrc = '88D VAD'
      CALL rdprof(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 3 .AND. suffix(1:3) == 'gps') THEN
 
      CALL rdgps(nx,ny,nz,nvar,nzua,mxua,nsrcua,uafname(ifile),        &
             anx,xs,ys,zp,su,sv,st,spres,shght,sqv,sqvsat,             &
             srcua,stnua,elevua,xua,yua,hgtua,obsua,                   &
             qualua,isrcua,nlevsua,                                    &
             rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 3 .AND. suffix(1:3) == 'air') THEN

      dfltsrc='AIRS'
      CALL rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 7 .AND. suffix(1:3) == 'AO_') THEN

      dfltsrc=suffix(1:7)
      CALL rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 7 .AND. suffix(1:3) == 'AL_') THEN

      dfltsrc=suffix(1:7)
      CALL rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    ELSE IF(lensuf == 7 .AND. suffix(1:3) == 'AM_') THEN

      dfltsrc=suffix(1:7)
      CALL rdraob(nvar,nzua,mxua,nsrcua,srcua,dfltsrc,uafname(ifile),  &
              stnua,elevua,xua,yua,hgtua,obsua,                        &
              qualua,isrcua,nlevsua,                                   &
              rmiss,nprev,ntotal,istatus)

    END IF
    nprev=ntotal
  END DO   ! ua file list
!
  nobsua=ntotal
!
!-----------------------------------------------------------------------
!
!  Set qobs based on source and height
!
!-----------------------------------------------------------------------
!
  DO ista=1,nobsua
    IF(isrcua(ista) > 0) THEN
      isrc=isrcua(ista)
      IF(nlevtab(isrc) < 2) THEN
        WRITE(6,'(a,a,i6,/a,i3)')                                       &
            ' Problem with error table for ',                           &
            ' upper-level data source',isrc,                            &
            '   nlevtab(isrc) =',nlevtab(isrc)
        CALL arpsstop("data problem",1)
      END IF
      DO ilev=1,nzua
        DO ktab=2,nlevtab(isrc)-1
          IF(hgtqsrc(ktab,isrc) > hgtua(ilev,ista)) EXIT
        END DO
!        126       CONTINUE
        wthi=  (hgtua(ilev,ista)-hgtqsrc(ktab-1,isrc))/                 &
             (hgtqsrc(ktab,isrc)-hgtqsrc(ktab-1,isrc))
        wthi=AMAX1(wthi,0.0)
        wthi=AMIN1(wthi,1.0)
        wtlo=1.0-wthi
        DO ivar=1,nvar
          qobsua(ivar,ilev,ista)=                                       &
                 wthi*qsrcua(ivar,ktab,  isrc) +                        &
                 wtlo*qsrcua(ivar,ktab-1,isrc)
        END DO
      END DO
    END IF
  END DO
  RETURN
END SUBROUTINE prepuaobs
