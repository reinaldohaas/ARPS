!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ANXITER                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE anxiter(nx,ny,nz,                                            &
           nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                     &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                         &
           mxpass,npass,iwstat,xs,ys,zs,icatg,xcor,nam_var,             &
           mpi_map,nmap,mpi_mapv,nmapv,                                 &
           xsng,ysng,hgtsng,thesng,trnsng,                              &
           obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,         &
           indexsng,np,ksng,ksngmax,indexua,kua,kuamax,                 &
           indexrad,krad,kradmax,                                       &
           xua,yua,hgtua,theua,trnua,                                   &
           obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,            &
           elvrad,xradc,yradc,                                          &
           distrad,uazmrad,vazmrad,hgtradc,theradc,trnradc,             &
           obsrad,odifrad,qobsrad,qualrad,                              &
           irad,isrcrad,nlevrad,ncolrad,ncolrad_mpi,                    &
           xretc,yretc,hgtretc,theretc,trnretc,                         &
           obsret,odifret,qobsret,qualret,                              &
           iret,isrcret,nlevret,ncolret,                                &
           srcsng,srcua,srcrad,srcret,                                  &
           ianxtyp,iusesng,iuseua,iuserad,iuseret,                      &
           xyrange,kpvrsq,wlim,zrange,zwlim,thrng,                      &
           trnropt,trnrcst,trnrng,rngsqi,knt,wgtsum,zsum,               &
           corsng,corua,corrad,corret,                                  &
           oanxsng,oanxua,oanxrad,oanxret,                              &
           anx,tem1,tem2,tem3,ibegin,iend,                              &
           tems2drsng,temr2drsng,tems3drua,temr3drua,                   &
           tems3dirad,temr3dirad,tems3drrad,temr3drrad,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Control ARPS analysis successive corrections iteration.
!  For each pass it calls routines to analyze at grid points,
!  analyze at obs points, computes new obs differences and
!  reports rms statistics.
!
!  On input the analyses at grid points and obs points
!  (anx, oanxsng and oanxua) should be initialized.
!  Either as as zero or with a first guess field and
!  that guess field interpolated to obs sites.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Version for OLAPS surface Keith Brewster, CAPS
!  March, 1994
!
!  3-D ARPS version Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!  Update for radar data, misc. improvements
!  Jan, 1996 KB
!
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!  July, 1996 KB
!
!  Changes made to anxiter argument list and calls to Bratseth routines
!  to accomodate changes made to speed up Bratseth routines.
!  Jan 5, 1997 KB
!
!  04/10/2003 (K. Brewster)
!  Updated print statements to remove format numbers and add flush
!  at the end of each iteration.
!
!  11/13/2003 (K. Brewster)
!  Merged with version that accomodates reduction in correlation due to
!  terrain variations for data near the surface.
!
!  2005 (K. W. Thomas)
!  Add MPI support.
!  Fix many bugs (mostly typos) that showed up when MPI'ing "adas".
!
!  11/05/2007 (K. W. Thomas)
!  Add MPI support for radial velocity processing.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx         Number of grid points in the x-direction (east/west)
!    ny         Number of grid points in the y-direction (north/south)
!    nz         Number of grid points in the vertical
!    nvar       Number of analysis variables
!    nvarradin  Number of variables in the obsrad array
!    nvarrad    Number of variables in the odifrad array
!    nzua       Maximumnumber of levels in UA observations
!    nzrdr      Maximum number of levels in a radar column
!    nzret      Maximum number of levels in a retrieval column
!    mxsng      Maximum number of single level data
!    mxua       Maximum number of upper air observations
!    mxrad      Maximum number of radars
!    mxcolrad   Maximum number of radar columns
!    mxcolret   Maximum number of retrieval columns
!    mxpass     Maximum number of iterations
!    npass      Number of iterations
!    iwstat     Status indicator for writing statistics
!
!    xs         x location of scalar pts (m)
!    ys         y location of scalar pts (m)
!    zs         z location of scalar pts (m)
!
!    mpi_map    Map of MPI communications
!    nmap       Number of entries in "mpi_map"
!    mpi_mapv   Map of MPI communications for radar radial velocity
!    nmapv      Number of entries in "mpi_mapv"
!
!    xsng       x location of single-level data
!    ysng       y location of single-level data
!    hgtsng     elevation of single-level data
!    thesng     theta (potential temperature) of single-level data
!    trnsng     model terrain at location of single-level data
!
!    obsng      single-level observations
!    odifsng    difference between single-level obs and analysis
!    qobsng     normalized observation error
!    qualsng    single-level data quality indicator
!    isrcsng    index of single-level data source
!    nobsng     number of single-level observations
!    indexsng   which processor "owns" which obs (MPI only)
!    np         number of processors (MPI only)
!    ksng       number of obs owned by each processor (MPI only)
!    ksngmax    largest "ksng" value (MPI only)
!
!    xua        x location of upper air data
!    yua        y location of upper air data
!    kua        number of multi-level obs owned by each processor (MPI only)
!    kuamax     largest "kua" value (MPI only)
!    hgtua      elevation of upper air data
!    theua      theta (potential temperature) of upper air data
!    trnua      model terrain at location of upper air data
!
!    obsua      upper air observations
!    odifua     difference between upper air obs and analysis
!    qobsua     normalized observation error
!    qualua     upper air data quality indicator
!    isrcua     index of upper air data source
!    nlevsua    number of levels of data for each upper air location
!    nobsua     number of upper air observations
!    indexua    which processor "owns" which obs (MPI only)
!
!    krad       number of multi-level obs owned by each processor (MPI only)
!    kradmax    largest "kua" value (MPI only)
!    indexrad   which processor "owns" which obs (MPI only)
!
!    anx        Analyzed variables (or first guess)
!    qback      Background (first guess) error
!
!    nradfil    number of radar files
!    fradname   file name for radar dataset
!    srcrad     name of radar sources
!
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!    xradc    x location of radar column
!    yradc    y location of radar column
!    irad     radar number
!    nlevrad  number of levels of radar data in each column
!    distrad  distance of radar column from source radar
!    uazmrad  u-component of radar beam for each column
!    vazmrad  v-component of radar beam for each column
!    hgtradc  height (m MSL) of radar observations
!    trnradc  model terrain at location of radar columns
!    obsrad   radar observations
!    oanxrad  analysis (first guess) value at radar data location
!    odifrad  difference between radar observation and analysis
!    qobsrad  normalized observation error
!    qualrad  radar data quality indicator
!    ncolrad  number of radar columns read-in, local if MPI
!    ncolrad_mpi  number of radar columns read-in, MPI only
!    istatus  status indicator
!
!    latret   latitude of retrieval radar  (degrees N)
!    lonret   longitude of retrieval radar (degrees E)
!    xretc    x location of retrieval column
!    yretc    y location of retrieval column
!    iret     retrieval number
!    nlevret  number of levels of retrieval data in each column
!    hgtretc  height (m MSL) of retrieval observations
!    trnretc  model terrain at location of retrieval observations
!    obsret   retrieval observations
!    odifret  difference between retrieval observation and analysis
!    qobsret  normalized observation error
!    qualret  retrieval data quality indicator
!    ncolret  number of retr columns read-in
!    istatus  status indicator
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
!-----------------------------------------------------------------------
!
!  Input Sizing Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nvarradin,nvarrad
  INTEGER :: nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  INTEGER :: mxpass,npass
!
!-----------------------------------------------------------------------
!
!  Input Grid Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  INTEGER :: icatg(nx,ny)
  REAL :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Input Observation Arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: nam_var(nvar)
  INTEGER :: nmap, nmapv
  INTEGER :: mpi_map(nmap,2), mpi_mapv(nmapv,2)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: thesng(mxsng)
  REAL :: trnsng(mxsng)
  REAL :: obsng(nvar,mxsng)
  REAL :: odifsng(nvar,mxsng)
  REAL :: qobsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
  INTEGER :: indexsng(mxsng)
  INTEGER :: np
  INTEGER :: ksng(np)
  INTEGER :: ksngmax
!
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  INTEGER :: kua(np)
  INTEGER :: kuamax
  REAL :: hgtua(nzua,mxua)
  REAL :: theua(nzua,mxua)
  REAL :: trnua(mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
  INTEGER :: indexua(mxua)
!
  REAL :: elvrad(mxrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: distrad(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: theradc(nzrdr,mxcolrad)
  REAL :: trnradc(mxcolrad)
  REAL :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: qobsrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
  INTEGER :: ncolrad_mpi
!
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: theretc(nzret,mxcolret)
  REAL :: trnretc(nzret,mxcolret)
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: qobsret(nvar,nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  MPI variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: krad(np)
  INTEGER :: kradmax
  INTEGER :: indexrad(mxcolrad)
!
!-----------------------------------------------------------------------
!
!  Input Analysis Control Variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcsng(nsrcsng)
  CHARACTER (LEN=8) :: srcua (nsrcua )
  CHARACTER (LEN=8) :: srcrad(nsrcrad)
  CHARACTER (LEN=8) :: srcret(nsrcret)

  INTEGER :: ianxtyp(mxpass)
  INTEGER :: iusesng(0:nsrcsng,mxpass)
  INTEGER :: iuseua(0:nsrcua,mxpass)
  INTEGER :: iuserad(0:nsrcrad,mxpass)
  INTEGER :: iuseret(0:nsrcret,mxpass)

  REAL :: xyrange(mxpass)
  REAL :: kpvrsq(nvar)
  REAL :: wlim
  REAL :: zrange(mxpass)
  REAL :: zwlim
  REAL :: thrng(mxpass)
  INTEGER :: trnropt(mxpass)
  REAL :: trnrcst(mxpass)
  REAL :: trnrng(mxpass)
  INTEGER :: iwstat
!
!-----------------------------------------------------------------------
!
!  Scratch Space
!
!-----------------------------------------------------------------------
!
  REAL :: rngsqi(nvar)
  INTEGER :: knt(nvar,nz)
  REAL :: wgtsum(nvar,nz)
  REAL :: zsum(nvar,nz)
!
!-----------------------------------------------------------------------
!
!  Output Variables at Observation Locations
!
!-----------------------------------------------------------------------
!
  REAL :: corsng(nvar,mxsng)
  REAL :: corua(nvar,nzua,mxua)
  REAL :: corrad(nvarrad,nzrdr,mxcolrad)
  REAL :: corret(nvar,nzret,mxcolret)
  REAL :: oanxsng(nvar,mxsng)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxret(nvar,nzret,mxcolret)
!
!-----------------------------------------------------------------------
!
!  Output Grid
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
!
!-----------------------------------------------------------------------
!
!  Work arrays, including MPI
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  INTEGER :: ibegin(ny)
  INTEGER :: iend(ny)

  REAL :: tems2drsng(nvar,ksngmax)
  REAL :: temr2drsng(nvar,ksngmax)
  REAL :: tems3drua(nvar,nzua,kuamax)
  REAL :: temr3drua(nvar,nzua,kuamax)
  REAL :: tems3dirad(nvarrad,nzrdr,kradmax)
  REAL :: temr3dirad(nvarrad,nzrdr,kradmax)
  REAL :: tems3drrad(nvarrad,nzrdr,kradmax)
  REAL :: temr3drrad(nvarrad,nzrdr,kradmax)
!
!-----------------------------------------------------------------------
!
!  Return status
!
!-----------------------------------------------------------------------
!
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL :: ftabinv,setexp
  INTEGER :: i,j,k,ipass,isrc,sngsw,uasw,radsw,retsw
  REAL :: rpass,zrngsq,thrngsq,trnrngsq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Note in the following the analysis scratch arrays
!  are used to hold the output statistics
!
!-----------------------------------------------------------------------
!
  ipass=0
  IF (myproc == 0) WRITE(6,'(a,i4)') 'Computing orig difference stats'
  CALL diffnst3d(nvar,                                                  &
                 nzua,nzret,mxsng,mxua,mxcolret,                        &
                 obsng,odifsng,oanxsng,isrcsng,qualsng,nobsng,          &
                 obsua,odifua,oanxua,isrcua,qualua,nlevsua,nobsua,      &
                 obsret,odifret,oanxret,iret,                           &
                 qualret,nlevret,ncolret,                               &
                 indexsng,indexua,                                      &
                 knt,wgtsum,zsum)
  IF (mp_opt > 0) THEN
    IF (nobsng > 0)                                                     &
      CALL mpi_2dr_collect(odifsng, nvar, mxsng, nobsng, indexsng, np,  &
        ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
    IF (nobsua > 0)                                                     &
      CALL mpi_3dr_collect(odifua, nvar, nzua, mxua, nobsua, indexua,   &
        np, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
  END IF

  IF (myproc == 0) THEN                                            
    WRITE(6,'(/a,i4,/a)')'  Statistics for analysis pass ',ipass,       &
                 '  ivar  name     knt    bias        rms'

    DO k=1,nvar
      WRITE(6,'(i6,2X,a6,i6,g12.3,g11.3)')                              &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
    END DO

    IF(iwstat > 0) THEN
      WRITE(iwstat,'(/a,i4,/a)')'  Statistics for analysis pass ',ipass,&
                 '  ivar  name     knt    bias        rms'
      DO k=1,nvar
        WRITE(iwstat,'(i6,2X,a6,i6,g12.3,g11.3)')                       &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
      END DO
    END IF
  END IF

  CALL diffnstra(nvar,nvarrad,nvarradin,nzrdr,mxrad,mxcolrad,           &
            elvrad,distrad,hgtradc,                                     &
            uazmrad,vazmrad,qualrad,irad,                               &
            obsrad,odifrad,oanxrad,                                     &
            nlevrad,ncolrad,indexrad,                                   &
            knt,wgtsum,zsum)

  IF (mp_opt > 0) THEN
    IF (ncolrad_mpi > 0) THEN
      CALL mpi_3dr_collect(odifrad, nvarrad, nzrdr, mxcolrad, ncolrad,  &
        indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3drrad, temr3drrad)
      CALL mpi_3di_collect(qualrad, nvarrad, nzrdr, mxcolrad, ncolrad,  &
        indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3dirad, temr3dirad)
    END IF
  END IF

  IF (myproc == 0) THEN
    WRITE(6,'(/a,i4,/a)') '  Radar data variables for pass ',ipass,     &
                 '  ivar  name    knt      bias        rms'
    DO k=1,nvarrad
      WRITE(6,'(i6,2X,a6,i7,g12.3,g11.3)')                              &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
    END DO

    IF(iwstat > 0) THEN
      WRITE(iwstat,'(/a,i4,/a)') '  Radar data variables for pass ',ipass,&
                 '  ivar  name    knt      bias        rms'
      DO k=1,nvarrad
        WRITE(iwstat,'(i6,2X,a6,i7,g12.3,g11.3)')                       &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
      END DO
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize correlation lookup table.
!
!-----------------------------------------------------------------------
!
  ftabinv=setexp(16.)
!
!-----------------------------------------------------------------------
!
!  Loop through npass analysis iterations
!
!-----------------------------------------------------------------------
!
  DO ipass=1,npass
!
!-----------------------------------------------------------------------
!
!  Set single-level usage switch based in iusesng
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) THEN
      WRITE(6,'(/a,i4)') ' Source usage switches for pass ',ipass
      WRITE(6,'(/a,i4)') '   Single level data'
    END IF
    sngsw=0
    DO isrc=1,nsrcsng
      IF(iusesng(isrc,ipass) > 0) THEN
        sngsw=1
        IF (myproc == 0 .AND. srcsng(isrc) /= 'NULL')                  &
           WRITE(6,'(a,a)') '      Using     ',srcsng(isrc)
      END IF
    END DO
    IF(myproc == 0 .AND. sngsw == 0) WRITE(6,'(a)') '      none'
!
!-----------------------------------------------------------------------
!
!  Set multiple-level usage switch based in iuseua
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(/a,i4)') '   Multiple level data'
    uasw=0
    DO isrc=1,nsrcua
      IF(iuseua(isrc,ipass) > 0) THEN
        uasw=1
        IF (myproc == 0 .AND. srcua(isrc) /= 'NULL')                   &
           WRITE(6,'(a,a)') '      Using     ',srcua(isrc)
      END IF
    END DO
    IF(myproc == 0 .AND. uasw == 0) WRITE(6,'(a)') '      none'
!
!-----------------------------------------------------------------------
!
!  Set radar-data usage switch based in iuserad
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(/a,i4)') '   Raw radar data'
    radsw=0
    DO isrc=1,nsrcrad
      IF(iuserad(isrc,ipass) > 0) THEN
        radsw=1
        IF (myproc == 0 .AND. srcrad(isrc) /= 'NULL')                  &
           WRITE(6,'(a,a)') '      Using     ',srcrad(isrc)
      END IF
    END DO
    IF(myproc == 0 .AND. radsw == 0) WRITE(6,'(a)') '      none'
!
!-----------------------------------------------------------------------
!
!  Set ret usage switch based in iuseret
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0) WRITE(6,'(/a,i4)') '   Retrieval radar data'
    retsw=0
    DO isrc=1,nsrcret
      IF(iuseret(isrc,ipass) > 0) THEN
        retsw=1
        IF (mp_opt > 0) THEN
          IF (myproc == 0) write(6,'(/a)') ' Retrival code is not MPI ready.'
          CALL arpsstop("MPI and retrieval not possible now",1)
        END IF
        IF (myproc == 0 .AND. srcret(isrc) /= 'NULL')                  &
           WRITE(6,'(a,a)') '      Using     ',srcret(isrc)
      END IF
    END DO
    IF(myproc == 0 .AND. retsw == 0) WRITE(6,'(a)') '      none'

    rpass=xyrange(ipass)*xyrange(ipass)
    zrngsq=zrange(ipass)*zrange(ipass)
    thrngsq=thrng(ipass)*thrng(ipass)
    trnrngsq=trnrng(ipass)*trnrng(ipass)
!
    IF (myproc == 0) THEN
      write(6,'(a,i3)') ' trnropt = ',trnropt(ipass)
      write(6,'(a,f10.1,a)') ' trnrcst = ',trnrcst(ipass),' m'
      write(6,'(a,f10.3)') ' trnrng = ',trnrng(ipass)
    END IF
!
    IF(ianxtyp(ipass) == 11) THEN

      IF (myproc == 0)                                                  &
        WRITE(6,'(/a,a,i4)') 'Barnes on grid using height',             &
                        ' ipass=',ipass
      CALL bargrd3d(nx,ny,nz,                                           &
             nvar,nvarrad,nzua,nzrdr,nzret,                             &
             mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                  &
             nsrcsng,nsrcua,nsrcrad,nsrcret,xs,ys,zs,                   &
             odifsng,xsng,ysng,hgtsng,qualsng,isrcsng,nobsng,           &
             odifua,xua,yua,hgtua,qualua,isrcua,nlevsua,nobsua,         &
             odifrad,xradc,yradc,hgtradc,                               &
             uazmrad,vazmrad,qualrad,                                   &
             irad,isrcrad,nlevrad,ncolrad,                              &
             odifret,xretc,yretc,hgtretc,                               &
             qualret,iret,isrcret,nlevret,ncolret,                      &
             iusesng(0,ipass),iuseua(0,ipass),                          &
             iuserad(0,ipass),iuseret(0,ipass),                         &
             sngsw,uasw,radsw,retsw,kpvrsq,rpass,                       &
             wlim,zrngsq,zwlim,                                         &
             rngsqi,knt,wgtsum,zsum,                                    &
             anx,istatus)
!
      IF (myproc == 0)                                                  &
        WRITE(6,'(/a,a,i4)') 'Barnes at obs sites using height',        &
                        ' ipass=',ipass
      CALL barobs3d(nvar,nvarrad,nzua,nzrdr,nzret,                      &
             mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                  &
             nsrcsng,nsrcua,nsrcrad,nsrcret,                            &
             odifsng,xsng,ysng,hgtsng,qualsng,isrcsng,nobsng,           &
             odifua,xua,yua,hgtua,qualua,isrcua,nlevsua,nobsua,         &
             odifrad,xradc,yradc,hgtradc,uazmrad,vazmrad,               &
             qualrad,irad,isrcrad,nlevrad,ncolrad,                      &
             odifret,xretc,yretc,hgtretc,                               &
             qualret,iret,isrcret,nlevret,ncolret,                      &
             iusesng(0,ipass),iuseua(0,ipass),                          &
             iuserad(0,ipass),iuseret(0,ipass),                         &
             sngsw,uasw,radsw,retsw,kpvrsq,rpass,                       &
             wlim,zrngsq,zwlim,                                         &
             rngsqi,knt,wgtsum,zsum,                                    &
             oanxsng,oanxua,oanxrad,oanxret,istatus)
!
    ELSE IF(ianxtyp(ipass) == 21) THEN
!
      CALL obscor(nvar,nvarrad,nzua,nzrdr,nzret,                        &
             mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                  &
             nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                  &
             qobsng,xsng,ysng,hgtsng,trnsng,                            &
             qualsng,isrcsng,icatsng,nobsng,                            &
             qobsua,xua,yua,hgtua,trnua,                                &
             qualua,isrcua,nlevsua,nobsua,                              &
             qobsrad,xradc,yradc,hgtradc,trnradc,                       &
             uazmrad,vazmrad,                                           &
             qualrad,irad,isrcrad,nlevrad,ncolrad,                      &
             qobsret,xretc,yretc,hgtretc,trnretc,                       &
             qualret,iret,isrcret,nlevret,ncolret,                      &
             iusesng(0,ipass),iuseua(0,ipass),                          &
             iuserad(0,ipass),iuseret(0,ipass),                         &
             sngsw,uasw,radsw,retsw,kpvrsq,rpass,rngsqi,                &
             wlim,zrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,        &
             corsng,corua,corrad,corret,istatus)

      IF (mp_opt > 0) THEN
        IF (nobsng > 0)                                                  &
          CALL mpi_2dr_collect(corsng, nvar, mxsng, nobsng, indexsng, np,&
            ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
        IF (nobsua > 0)                                                  &
          CALL mpi_3dr_collect(corua, nvar, nzua, mxua, nobsua, indexua, &
            np, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
        IF (ncolrad_mpi > 0)                                             &
          CALL mpi_3dr_collect(corrad, nvarrad, nzrdr, mxcolrad, ncolrad,&
            indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3drrad,    &
            temr3drrad)
      END IF

!
!-----------------------------------------------------------------------
!
!  Calculate covariance sums at observation points for
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0)                                                  &
        WRITE(6,'(/a,a,i4)') 'Bratseth on grid using height',           &
                        ' ipass=',ipass
      CALL brtgrd3d(nx,ny,nz,nvar,nvarrad,                              &
            nzua,nzrdr,nzret,                                           &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xs,ys,zs,icatg,xcor,    &
            odifsng,corsng,xsng,ysng,hgtsng,trnsng,                     &
            isrcsng,icatsng,nobsng,                                     &
            odifua,corua,xua,yua,hgtua,trnua,isrcua,nlevsua,nobsua,     &
            odifrad,corrad,xradc,yradc,hgtradc,trnradc,                 &
            uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,               &
            odifret,corret,xretc,yretc,hgtretc,trnretc,                 &
            iret,isrcret,nlevret,ncolret,                               &
            iusesng(0,ipass),iuseua(0,ipass),                           &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            sngsw,uasw,radsw,retsw,kpvrsq,rpass,                        &
            wlim,zrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,rngsqi,  &
            anx,tem1,tem2,tem3,ibegin,iend,istatus)
!
      IF (myproc == 0)                                                  &
        WRITE(6,'(a,a,i4)') 'Bratseth at obs sites using height',       &
                        ' ipass=',ipass
      CALL brtobs3d(nvar,nvarrad,nzua,nzrdr,nzret,                      &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                   &
            odifsng,qobsng,corsng,xsng,ysng,hgtsng,trnsng,              &
            isrcsng,icatsng,nobsng,                                     &
            odifua,qobsua,corua,xua,yua,hgtua,trnua,                    &
            isrcua,nlevsua,nobsua,                                      &
            odifrad,qobsrad,corrad,xradc,yradc,hgtradc,trnradc,         &
            uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,               &
            odifret,qobsret,corret,xretc,yretc,hgtretc,trnretc,         &
            iret,isrcret,nlevret,ncolret,                               &
            iusesng(0,ipass),iuseua(0,ipass),                           &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            sngsw,uasw,radsw,retsw,kpvrsq,rpass,                        &
            wlim,zrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,rngsqi,  &
            oanxsng,oanxua,oanxrad,oanxret,istatus)

    ELSE IF(ianxtyp(ipass) == 12) THEN

      IF (myproc == 0)                                                  &
        WRITE(6,'(/a)') 'Method ianxtyp = 12 is still not implemented.'
      CALL arpsstop('Not implemented method called inside anxiter.',1)
      stop

    ELSE   ! ianxtyp(ipass) == 22
!
      CALL obscor(nvar,nvarrad,nzua,nzrdr,nzret,                        &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                   &
            qobsng,xsng,ysng,thesng,trnsng,                             &
            qualsng,isrcsng,icatsng,nobsng,                             &
            qobsua,xua,yua,theua,trnua,                                 &
            qualua,isrcua,nlevsua,nobsua,                               &
            qobsrad,xradc,yradc,theradc,trnradc,                        &
            uazmrad,vazmrad,                                            &
            qualrad,irad,isrcrad,nlevrad,ncolrad,                       &
            qobsret,xretc,yretc,theretc,trnretc,                        &
            qualret,iret,isrcret,nlevret,ncolret,                       &
            iusesng(0,ipass),iuseua(0,ipass),                           &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            sngsw,uasw,radsw,retsw,kpvrsq,rpass,rngsqi,                 &
            wlim,thrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,        &
            corsng,corua,corrad,corret,istatus)

      IF (mp_opt > 0) THEN
        IF (nobsng > 0)                                                 &
          CALL mpi_2dr_collect(corsng, nvar, mxsng, nobsng, indexsng, np,&
            ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
        IF (nobsua > 0)                                                 &
          CALL mpi_3dr_collect(corua, nvar, nzua, mxua, nobsua, indexua,&
            np, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
        IF (ncolrad_mpi > 0)                                             &
          CALL mpi_3dr_collect(corrad, nvarrad, nzrdr, mxcolrad, ncolrad,&
            indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3drrad,    &
            temr3drrad)
      END IF

!
      IF (myproc == 0)                                                  &
        WRITE(6,'(/a,a,i4)') 'Bratseth on grid using theta',            &
                        ' ipass=',ipass
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem3(i,j,k)=anx(i,j,k,4)
          END DO
        END DO
      END DO
!
      CALL brtgrd3d(nx,ny,nz,nvar,nvarrad,                              &
            nzua,nzrdr,nzret,                                           &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                        &
            xs,ys,tem3,icatg,xcor,                                      &
            odifsng,corsng,xsng,ysng,thesng,trnsng,                     &
            isrcsng,icatsng,nobsng,                                     &
            odifua,corua,xua,yua,theua,trnua,                           &
            isrcua,nlevsua,nobsua,                                      &
            odifrad,corrad,xradc,yradc,theradc,trnradc,                 &
            uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,               &
            odifret,corret,xretc,yretc,theretc,trnretc,                 &
            iret,isrcret,nlevret,ncolret,                               &
            iusesng(0,ipass),iuseua(0,ipass),                           &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            sngsw,uasw,radsw,retsw,kpvrsq,rpass,                        &
            wlim,thrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,rngsqi, &
            anx,tem1,tem2,tem3,ibegin,iend,istatus)
!
      IF (myproc == 0)                                                  &
        WRITE(6,'(a,a,i4)') 'Bratseth at obs sites using theta',        &
                        ' ipass=',ipass
!
!-----------------------------------------------------------------------
!
      CALL brtobs3d(nvar,nvarrad,nzua,nzrdr,nzret,                      &
             mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                  &
             nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                  &
             odifsng,qobsng,corsng,xsng,ysng,thesng,trnsng,             &
             isrcsng,icatsng,nobsng,                                    &
             odifua,qobsua,corua,xua,yua,theua,trnua,                   &
             isrcua,nlevsua,nobsua,                                     &
             odifrad,qobsrad,corrad,xradc,yradc,theradc,trnradc,        &
             uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,              &
             odifret,qobsret,corret,xretc,yretc,theretc,trnretc,        &
             iret,isrcret,nlevret,ncolret,                              &
             iusesng(0,ipass),iuseua(0,ipass),                          &
             iuserad(0,ipass),iuseret(0,ipass),                         &
             sngsw,uasw,radsw,retsw,kpvrsq,rpass,                       &
             wlim,zrngsq,trnropt(ipass),trnrcst(ipass),trnrngsq,rngsqi, &
             oanxsng,oanxua,oanxrad,oanxret,istatus)
    END IF
!
!-----------------------------------------------------------------------
!
!  Note in the following the analysis scratch arrays
!  are used to hold the output statistics
!
!-----------------------------------------------------------------------
!
    IF (myproc == 0)                                                    &
      WRITE(6,'(a,i4)') 'Computing new differences, ipass=',ipass
    CALL diffnst3d(nvar,                                                &
                   nzua,nzret,mxsng,mxua,mxcolret,                      &
                   obsng,odifsng,oanxsng,isrcsng,qualsng,nobsng,        &
                   obsua,odifua,oanxua,isrcua,qualua,nlevsua,nobsua,    &
                   obsret,odifret,oanxret,iret,                         &
                   qualret,nlevret,ncolret,                             &
                   indexsng,indexua,                                    &
                   knt,wgtsum,zsum)
 
    IF (mp_opt > 0) THEN
      IF (nobsng > 0)                                                   &
        CALL mpi_2dr_collect(odifsng, nvar, mxsng, nobsng, indexsng, np,&
          ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
      IF (nobsua > 0)                                                   &
        CALL mpi_3dr_collect(odifua, nvar, nzua, mxua, nobsua, indexua, &
          np, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
    END IF

    IF (myproc == 0) THEN
      WRITE(6,'(/a,i4,/a)')'  Statistics for analysis pass ',ipass,     &
                 '  ivar  name     knt    bias        rms'
      DO k=1,nvar
        WRITE(6,'(i6,2X,a6,i6,g12.3,g11.3)')                            &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
      END DO
      IF(iwstat > 0) THEN
        WRITE(iwstat,'(/a,i4,/a)')'  Statistics for analysis pass ',ipass,&
                 '  ivar  name     knt    bias        rms'
        DO k=1,nvar
          WRITE(iwstat,'(i6,2X,a6,i6,g12.3,g11.3)')                     &
              k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
        END DO
      END IF
    END IF
 
    CALL diffnstra(nvar,nvarrad,nvarradin,nzrdr,mxrad,mxcolrad,         &
              elvrad,distrad,hgtradc,                                   &
              uazmrad,vazmrad,qualrad,irad,                             &
              obsrad,odifrad,oanxrad,                                   &
              nlevrad,ncolrad,indexrad,                                 &
              knt,wgtsum,zsum)

    IF (mp_opt > 0) THEN
      IF (ncolrad_mpi > 0) THEN
        CALL mpi_3dr_collect(odifrad, nvarrad, nzrdr, mxcolrad, ncolrad,&
          indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3drrad,     &
          temr3drrad)
        CALL mpi_3di_collect(qualrad, nvarrad, nzrdr, mxcolrad, ncolrad,&
          indexrad, np, krad, kradmax, mpi_mapv, nmapv, tems3dirad, temr3dirad)
      END IF
    END IF

    IF (myproc == 0) THEN
      WRITE(6,'(/a,i4,/a)') '  Radar data variables for pass ',ipass,   &
                 '  ivar  name    knt      bias        rms'
      DO k=1,nvarrad
        WRITE(6,'(i6,2X,a6,i7,g12.3,g11.3)')                            &
            k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
      END DO

      IF(iwstat > 0) THEN
        WRITE(iwstat,'(/a,i4,/a)') '  Radar data variables for pass ',  &
                 ipass,'  ivar  name    knt      bias        rms'
        DO k=1,nvarrad
          WRITE(iwstat,'(i6,2X,a6,i7,g12.3,g11.3)')                     &
                k,nam_var(k),knt(k,1),wgtsum(k,1),zsum(k,1)
        END DO
      END IF
   END IF
!
! Because each analysis iteration can take a significant amount
! time if there are many data to process, flush the
! stdout buffer so the progress can be tracked.
!
    CALL flush(6)
    IF(iwstat > 0) CALL flush(iwstat)

  END DO
!
  istatus=0
  RETURN
END SUBROUTINE anxiter
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFFNST3d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE diffnst3d(nvar,                                              &
           nzua,nzret,mxsng,mxua,mxcolret,                              &
           obsng,odifsng,oanxsng,isrcsng,qualsng,nobsng,                &
           obsua,odifua,oanxua,isrcua,qualua,nlevsua,nobsua,            &
           obsret,odifret,oanxret,iret,                                 &
           qualret,nlevret,ncolret,                                     &
           indexsng,indexua,                                            &
           knt,dmean,drms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate new observation differences
!  and report statistics on the differences.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  March, 1994
!
!  MODIFICATION HISTORY:
!  July, 1995  (Keith Brewster)
!  3D ARPS version
!
!  May, 1996 (KB)
!  Added retrieval data.
!
!  Jan, 2006 (KWT)
!  Add MPI support for stats.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nvar
  INTEGER :: nzua,nzret,mxsng,mxua,mxcolret
!
  REAL :: obsng(nvar,mxsng)
  REAL :: odifsng(nvar,mxsng)
  REAL :: oanxsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: nobsng
!
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: oanxua(nvar,nzua,mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: nobsua
!
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: oanxret(nvar,nzret,mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: ncolret
  INTEGER :: indexsng(mxsng)
  INTEGER :: indexua(mxua)
!
  INTEGER :: knt(nvar)
  REAL :: dmean(nvar)
  REAL :: drms(nvar)
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ivar,jsta,klev
  INTEGER :: k
  REAL :: flknt
!
!-----------------------------------------------------------------------
!
!  Double precision is required here so that MPI and non-MPI runs can
!  produce the same answers.
!
!-----------------------------------------------------------------------
!
  DOUBLE PRECISION dmean2(5)

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO ivar=1,nvar
    knt(ivar)=0
    dmean(ivar)=0.
    dmean2(ivar)=0.
    drms(ivar)=0.
  END DO
!
!-----------------------------------------------------------------------
!
!  Find differences and form sums at single-level sites
!
!-----------------------------------------------------------------------
!
  DO jsta=1,nobsng
  IF (mp_opt == 0 .OR. (mp_opt > 0 .AND. indexsng(jsta) == myproc)) THEN
!
!-----------------------------------------------------------------------
!
!    We have to check "isrcsng" to make sure we're a valid station to
!    use.  Normally "isrc*" is 0 for no data and positive otherwise.
!    "isrcsng" is the exception.  "-1" is used to flag data that can
!    never be used by the cloud soundings (cloudopt=1).  Data with
!    "isrc" of zero is outside of the domain, but still could be used
!    for cloud soundings.  "-1" data is SUPROB joined, so originals
!    data is no longer used.
!
!-----------------------------------------------------------------------
!
      IF(isrcsng(jsta) > 0) THEN
        DO ivar=1,nvar
          IF(qualsng(ivar,jsta) > 0) THEN
            odifsng(ivar,jsta)=obsng(ivar,jsta)-oanxsng(ivar,jsta)
            knt(ivar)=knt(ivar)+1

            dmean(ivar)=dmean(ivar)-odifsng(ivar,jsta)
            dmean2(ivar)=dmean2(ivar)-odifsng(ivar,jsta)
            drms(ivar)=drms(ivar)+                                        &
                       odifsng(ivar,jsta)*odifsng(ivar,jsta)
          END IF
        END DO
      END IF
  END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Find differences at multiple-level sites
!  Add to sums
!
!-----------------------------------------------------------------------
!
  DO jsta=1,nobsua
  IF (mp_opt == 0 .OR. (mp_opt > 0 .AND. indexua(jsta) == myproc)) THEN
      IF(isrcua(jsta) > 0) THEN
        DO klev=1,nlevsua(jsta)
          DO ivar=1,nvar
            IF(qualua(ivar,klev,jsta) > 0) THEN
              odifua(ivar,klev,jsta)=                                     &
                  obsua(ivar,klev,jsta)-oanxua(ivar,klev,jsta)
              knt(ivar)=knt(ivar)+1
              dmean(ivar)=dmean(ivar)-odifua(ivar,klev,jsta)
              dmean2(ivar)=dmean2(ivar)-odifua(ivar,klev,jsta)
              drms(ivar)=drms(ivar)+                                      &
                  odifua(ivar,klev,jsta)*odifua(ivar,klev,jsta)
            END IF
          END DO
        END DO
      END IF
  END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Find differences at retrieval columns
!  Add to sums
!
!-----------------------------------------------------------------------
!
  DO jsta=1,ncolret
    IF(iret(jsta) > 0) THEN
      DO klev=1,nlevret(jsta)
        DO ivar=1,nvar
          IF(qualret(ivar,klev,jsta) > 0) THEN
            odifret(ivar,klev,jsta)=                                    &
                obsret(ivar,klev,jsta)-oanxret(ivar,klev,jsta)
            knt(ivar)=knt(ivar)+1
            dmean(ivar)=dmean(ivar)-odifret(ivar,klev,jsta)
            dmean2(ivar)=dmean2(ivar)-odifret(ivar,klev,jsta)
            drms(ivar)=drms(ivar)+                                      &
                odifret(ivar,klev,jsta)*odifret(ivar,klev,jsta)
          END IF
        END DO
      END DO
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate stats from the sums formed above
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    DO k=1,nvar
      CALL mpsumi(knt(k),1)
      CALL mpsumr(dmean(k),1)
      CALL mpsumdp(dmean2(k),1)
      CALL mpsumr(drms(k),1)
    END DO
  END IF
  DO ivar=1,nvar
    IF(knt(ivar) > 0) THEN
      flknt=FLOAT(knt(ivar))
      dmean(ivar)=dmean(ivar)/flknt
      dmean2(ivar)=dmean2(ivar)/flknt
      dmean(ivar) = dmean2(ivar)
      drms(ivar)=SQRT(drms(ivar)/flknt)
    ELSE
      dmean(ivar)=0.
      drms(ivar)=0.
    END IF
  END DO
!
  RETURN
END SUBROUTINE diffnst3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFFNSTRA                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE diffnstra(nvar,nvarrad,nvarradin,nzrdr,mxrad,mxcolrad,       &
           elvrad,distrad,hgtradc,                                      &
           uazmrad,vazmrad,qualrad,irad,                                &
           obsrad,odifrad,oanxrad,                                      &
           nlevrad,ncolrad,indexrad,                                    &
           knt,dmean,drms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate radar observation differences
!  and report statistics on the differences.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  September, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (KB)
!  Added radar data and other improvements.
!
!  August, 2001 (KB)
!  Corrected call to dhdr which now returns dhdr not local
!  elevation angle.  Modified calculation of odifrad to speed
!  convergence to vr.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nvar,nvarrad,nvarradin,nzrdr,mxrad,mxcolrad
!
  REAL :: elvrad(mxrad)
  REAL :: distrad(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: irad(mxcolrad)
  REAL :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: ncolrad
  INTEGER :: indexrad(mxcolrad)
!
  INTEGER :: knt(nvar)
  REAL :: dmean(nvar)
  REAL :: drms(nvar)
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: icol,ivar,kr,jrad
  INTEGER :: k
  REAL :: pi,dtr,delz,eleva,range,flknt
  REAL :: vr,vrdiff,avrdif,dhdr,dsdr,dsdrinv
  REAL :: vr_miss
!
!-----------------------------------------------------------------------
!
!  QC and conversion parameters
!
!  vrqcthr   maximum radial velocity diff from background to use
!  vrcnvthr  lower limit of the dot product between the radial
!            direction and the analysis wind component to use data
!
!-----------------------------------------------------------------------
!
  REAL :: vrqcthr,vrcnvthr
  PARAMETER (vrqcthr=20.,   & ! maximum vrdiff from background to use
         vrcnvthr=-0.2) ! lower limit of the dot product to use data
!
!-----------------------------------------------------------------------
!
!  Double precision is required here so that MPI and non-MPI runs can
!  produce the same answers.
!
!-----------------------------------------------------------------------
!
  DOUBLE PRECISION dmean2(5)

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Some initializations
!
!-----------------------------------------------------------------------
!
  vr_miss=-999.
  pi=4.*ATAN(1.0)
  dtr=pi/180.
!
  DO ivar=1,nvar
    knt(ivar)=0
    dmean(ivar)=0.
    dmean2(ivar)=0.
    drms(ivar)=0.
  END DO
!
!-----------------------------------------------------------------------
!
!  Find radial velocity differences at radar data points
!  Add to sums
!
!-----------------------------------------------------------------------
!
  DO icol=1,ncolrad
   IF (mp_opt == 0 .OR. (mp_opt > 0 .AND. indexrad(icol) == myproc)) THEN
    IF(irad(icol) > 0) THEN
      jrad=irad(icol)
!
      DO kr=1,nlevrad(icol)
        IF(qualrad(1,kr,icol) > 0 .AND. qualrad(2,kr,icol) > 0 ) THEN
!
          delz=hgtradc(kr,icol)-elvrad(jrad)
          CALL beamelv(delz,distrad(icol),eleva,range)
          CALL dhdrange(eleva,range,dhdr)
          dsdr=SQRT(MAX(0.,(1.-dhdr*dhdr)))
          IF(dsdr /= 0.) THEN
            dsdrinv=1./dsdr
          ELSE
            dsdrinv=0.
          END IF
!
          vr=(uazmrad(icol)*oanxrad(1,kr,icol) +                        &
              vazmrad(icol)*oanxrad(2,kr,icol)) * dsdr
          vrdiff=obsrad(2,kr,icol)-vr
          avrdif=ABS(vrdiff)

          IF(avrdif < vrqcthr) THEN
!
!-----------------------------------------------------------------------
!
!  Find increment in analysis u and v needed to make vrdiff zero
!
!-----------------------------------------------------------------------
!
            odifrad(1,kr,icol)=dsdrinv*uazmrad(icol)*vrdiff
            odifrad(2,kr,icol)=dsdrinv*vazmrad(icol)*vrdiff
!
            knt(1)=knt(1)+1
            dmean(1)=dmean(1)-odifrad(1,kr,icol)
            dmean2(1)=dmean2(1)-odifrad(1,kr,icol)
            drms(1)=drms(1)+                                            &
                odifrad(1,kr,icol)*odifrad(1,kr,icol)

            knt(2)=knt(2)+1
            dmean(2)=dmean(2)-odifrad(2,kr,icol)
            dmean2(2)=dmean2(2)-odifrad(2,kr,icol)
            drms(2)=drms(2)+                                            &
                odifrad(2,kr,icol)*odifrad(2,kr,icol)

          ELSE

            odifrad(1,kr,icol)=vr_miss
            qualrad(1,kr,icol)=-199
            odifrad(2,kr,icol)=vr_miss
            qualrad(2,kr,icol)=-199

          END IF
        END IF
        IF(qualrad(3,kr,icol) > 0) THEN
          odifrad(3,kr,icol)=obsrad(1,kr,icol)-oanxrad(3,kr,icol)
!          print *, ' kr,icol,obsrad,oanx= ',kr,icol,
!    +                 obsrad(1,kr,icol),oanxrad(3,kr,icol)
          knt(5)=knt(5)+1
          dmean(5)=dmean(5)-odifrad(3,kr,icol)
          dmean2(5)=dmean2(5)-odifrad(3,kr,icol)
          drms(5)=drms(5)+                                              &
              odifrad(3,kr,icol)*odifrad(3,kr,icol)
        END IF
      END DO
    END IF
   END IF
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate stats from the sums formed above
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    DO k=1,nvar
      CALL mpsumi(knt(k),1)
      CALL mpsumr(dmean(k),1)
      CALL mpsumdp(dmean2(k),1)
      CALL mpsumr(drms(k),1)
    END DO
  END IF
  DO ivar=1,nvar
    IF(knt(ivar) > 0) THEN
      flknt=FLOAT(knt(ivar))
      dmean(ivar)=dmean(ivar)/flknt
      dmean2(ivar)=dmean2(ivar)/flknt
      dmean(ivar) = dmean2(ivar)
      drms(ivar)=SQRT(drms(ivar)/flknt)
    ELSE
      dmean(ivar)=0.
      drms(ivar)=0.
    END IF
  END DO
!
  RETURN
END SUBROUTINE diffnstra
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BAROBS3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE barobs3d(nvar,nvarrad,nzua,nzrdr,nzret,                      &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,                              &
           odifsng,xsng,ysng,hgtsng,qualsng,isrcsng,nobsng,             &
           odifua,xua,yua,hgtua,qualua,isrcua,nlevsua,nobsua,           &
           odifrad,xradc,yradc,hgtradc,uazmrad,vazmrad,                 &
           qualrad,irad,isrcrad,nlevrad,ncolrad,                        &
           odifret,xretc,yretc,hgtretc,                                 &
           qualret,iret,isrcret,nlevret,ncolret,                        &
           iusesng,iuseua,iuserad,iuseret,                              &
           sngsw,uasw,radsw,retsw,kpvrsq,rpass,                         &
           wlim,zrngsq,zwlim,                                           &
           rngsqi,knt,wgtsum,zsum,                                      &
           oanxsng,oanxua,oanxrad,oanxret,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Objectively analyze data at observation pts
!  Barnes weight function is employed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (KB)
!  Added radar data and other improvements.
!
!  July, 1996 (KB)
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Sizing arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nvar,nvarrad,nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: odifsng(nvar,mxsng),                                          &
       xsng(mxsng),                                                     &
       ysng(mxsng),                                                     &
       hgtsng(mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: nobsng
  REAL :: odifua(nvar,nzua,mxua)             ! variable to analyse
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Scratch space
!
!-----------------------------------------------------------------------
!
  INTEGER :: knt(nvar)
  REAL :: rngsqi(nvar)
  REAL :: wgtsum(nvar)
  REAL :: zsum(nvar)
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)
  INTEGER :: sngsw,uasw,radsw,retsw
  REAL :: kpvrsq(nvar)
  REAL :: rpass
  REAL :: wlim
  REAL :: zrngsq
  REAL :: zwlim
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  REAL :: oanxsng(nvar,mxsng)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxret(nvar,nzret,mxcolret)
  INTEGER :: istatus
include 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL :: texp
  REAL :: rlimsq,zlimsq,zsqinv,trnsqinv
  REAL :: dist,dzsq,wgt
  INTEGER :: ista,jsta,k,klev,ivar
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rlimsq=-rpass*ALOG(wlim)
  DO ivar=1,nvar
    rngsqi(ivar)=1./(kpvrsq(ivar)*rpass)
  END DO
!
  zlimsq=-zrngsq*ALOG(zwlim)
  zsqinv=1./zrngsq
!
!-----------------------------------------------------------------------
!
!  Loop over single stations.
!
!-----------------------------------------------------------------------
!
  DO ista=1,nobsng
   IF(isrcsng(ista) > 0) THEN
    IF(iusesng(isrcsng(ista)) > 0) THEN
      DO ivar=1,nvar
        knt(ivar)=0
        wgtsum(ivar)=0.
        zsum(ivar)=0.
      END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
      IF(sngsw > 0) THEN
        DO jsta=1,nobsng
         IF(isrcsng(jsta) > 0) THEN
          IF(iusesng(isrcsng(jsta)) > 0) THEN
            dist=(xsng(ista)-xsng(jsta))*(xsng(ista)-xsng(jsta))        &
                +(ysng(ista)-ysng(jsta))*(ysng(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dzsq=(hgtsng(ista)-hgtsng(jsta)) *                        &
                   (hgtsng(ista)-hgtsng(jsta))
              IF(dzsq < zlimsq) THEN
                DO ivar=1,nvar
                  IF(qualsng(ivar,jsta) > 0) THEN
                    knt(ivar)=knt(ivar)+1
                    wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                    wgtsum(ivar)=wgtsum(ivar) + wgt
                    zsum(ivar)=zsum(ivar)+                              &
                          wgt*odifsng(ivar,jsta)
                  END IF
                END DO
              END IF
            END IF     ! horizontal dist check
          END IF     ! source check
         END IF      ! isrcsng check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
      IF(uasw > 0) THEN
        DO jsta=1,nobsua
          IF(iuseua(isrcua(jsta)) > 0) THEN
            dist=(xsng(ista)-xua(jsta))*(xsng(ista)-xua(jsta))          &
                +(ysng(ista)-yua(jsta))*(ysng(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevsua(jsta)
                dzsq=(hgtsng(ista)-hgtua(klev,jsta)) *                  &
                     (hgtsng(ista)-hgtua(klev,jsta))
                IF(dzsq < zlimsq) THEN
                  DO ivar=1,nvar
                    IF(qualua(ivar,klev,jsta) > 0) THEN
                      knt(ivar)=knt(ivar)+1
                      wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                      wgtsum(ivar)=wgtsum(ivar) + wgt
                      zsum(ivar)=zsum(ivar)+                            &
                          wgt*odifua(ivar,klev,jsta)
                    END IF
                  END DO
                END IF
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
      IF(radsw > 0) THEN
        DO jsta=1,ncolrad
          IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
            dist=(xsng(ista)-xradc(jsta))*(xsng(ista)-xradc(jsta))      &
                +(ysng(ista)-yradc(jsta))*(ysng(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevrad(jsta)
                dzsq=(hgtsng(ista)-hgtradc(klev,jsta)) *                &
                     (hgtsng(ista)-hgtradc(klev,jsta))
                IF(dzsq < zlimsq) THEN
                  IF(qualrad(1,klev,jsta) > 0) THEN
                    knt(1)=knt(1)+1
                    wgt=texp(-dist*rngsqi(1) -dzsq*zsqinv)
                    wgtsum(1)=wgtsum(1)+wgt*                            &
                              ABS(uazmrad(jsta))
                    zsum(1)=zsum(1)+                                    &
                         wgt*odifrad(1,klev,jsta) *                     &
                                ABS(uazmrad(jsta))
                  END IF
                  IF(qualrad(2,klev,jsta) > 0) THEN
                    knt(2)=knt(2)+1
                    wgt=texp(-dist*rngsqi(2) -dzsq*zsqinv)
                    wgtsum(2)=wgtsum(2)+wgt*                            &
                              ABS(vazmrad(jsta))
                    zsum(2)=zsum(2)+                                    &
                         wgt*odifrad(2,klev,jsta) *                     &
                                ABS(vazmrad(jsta))
                  END IF
                  IF(qualrad(3,klev,jsta) > 0) THEN
                    knt(5)=knt(5)+1
                    wgt=texp(-dist*rngsqi(5) -dzsq*zsqinv)
                    wgtsum(5)=wgtsum(5)+wgt
                    zsum(5)=zsum(5)+                                    &
                         wgt*odifrad(3,klev,jsta)
                  END IF
                END IF
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval obs
!
!-----------------------------------------------------------------------
!
      IF(retsw > 0) THEN
        DO jsta=1,ncolret
          IF(iuseret(isrcret(iret(jsta))) > 0) THEN
            dist=(xsng(ista)-xretc(jsta))*(xsng(ista)-xretc(jsta))      &
                +(ysng(ista)-yretc(jsta))*(ysng(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevret(jsta)
                dzsq=(hgtsng(ista)-hgtretc(klev,jsta)) *                &
                     (hgtsng(ista)-hgtretc(klev,jsta))
                IF(dzsq < zlimsq) THEN
                  DO ivar=1,nvar
                    IF(qualret(ivar,klev,jsta) > 0) THEN
                      knt(ivar)=knt(ivar)+1
                      wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                      wgtsum(ivar)=wgtsum(ivar) + wgt
                      zsum(ivar)=zsum(ivar)+                            &
                          wgt*odifret(ivar,klev,jsta)
                    END IF
                  END DO
                END IF
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Divide by sum of weights
!
!-----------------------------------------------------------------------
!
      DO ivar=1,nvar
        IF(wgtsum(ivar) > wlim) oanxsng(ivar,ista)=oanxsng(ivar,ista) + &
                           zsum(ivar)/wgtsum(ivar)
      END DO
    END IF
   END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop over all multi-level obs
!
!-----------------------------------------------------------------------
!
  DO ista=1,nobsua
    IF(iuseua(isrcua(ista)) > 0) THEN
      DO k=1,nlevsua(ista)
        DO ivar=1,nvar
          zsum(ivar)=0.
          wgtsum(ivar)=0.
          knt(ivar)=0
        END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
        IF(sngsw > 0) THEN
          DO jsta=1,nobsng
           IF(isrcsng(jsta) > 0) THEN
            IF(iusesng(isrcsng(jsta)) > 0) THEN
              dist=(xua(ista)-xsng(jsta))*(xua(ista)-xsng(jsta))        &
                  +(yua(ista)-ysng(jsta))*(yua(ista)-ysng(jsta))
              IF(dist < rlimsq) THEN
                dzsq=(hgtua(k,ista)-hgtsng(jsta)) *                     &
                     (hgtua(k,ista)-hgtsng(jsta))
                IF(dzsq < zlimsq) THEN
                  DO ivar=1,nvar
                    IF(qualsng(ivar,jsta) > 0) THEN
                      knt(ivar)=knt(ivar)+1
                      wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                      wgtsum(ivar)=wgtsum(ivar) + wgt
                      zsum(ivar)=zsum(ivar)+                            &
                           wgt*odifsng(ivar,jsta)
                    END IF
                  END DO
                END IF
              END IF     ! horizontal dist check
            END IF     ! source check
           END IF      ! isrcsng check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
        IF(uasw > 0) THEN
          DO jsta=1,nobsua
            IF(iuseua(isrcua(jsta)) > 0) THEN
              dist=(xua(ista)-xua(jsta))*(xua(ista)-xua(jsta))          &
                  +(yua(ista)-yua(jsta))*(yua(ista)-yua(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevsua(jsta)
                  dzsq=(hgtua(k,ista)-hgtua(klev,jsta)) *               &
                       (hgtua(k,ista)-hgtua(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualua(ivar,klev,jsta) > 0) THEN
                        knt(ivar)=knt(ivar)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar)=wgtsum(ivar) + wgt
                        zsum(ivar)=zsum(ivar)+                          &
                            wgt*odifua(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
        IF(radsw > 0) THEN
          DO jsta=1,ncolrad
            IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
              dist=(xua(ista)-xradc(jsta))*(xua(ista)-xradc(jsta))      &
                  +(yua(ista)-yradc(jsta))*(yua(ista)-yradc(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevrad(jsta)
                  dzsq=(hgtua(k,ista)-hgtradc(klev,jsta)) *             &
                       (hgtua(k,ista)-hgtradc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    IF(qualrad(1,klev,jsta) > 0) THEN
                      knt(1)=knt(1)+1
                      wgt=texp(-dist*rngsqi(1) -dzsq*zsqinv)
                      wgtsum(1)=wgtsum(1)+wgt*                          &
                               ABS(uazmrad(jsta))
                      zsum(1)=zsum(1)+                                  &
                           wgt*odifrad(1,klev,jsta) *                   &
                               ABS(uazmrad(jsta))
                    END IF
                    IF(qualrad(2,klev,jsta) > 0) THEN
                      knt(2)=knt(2)+1
                      wgt=texp(-dist*rngsqi(2) -dzsq*zsqinv)
                      wgtsum(2)=wgtsum(2)+wgt*                          &
                               ABS(vazmrad(jsta))
                      zsum(2)=zsum(2)+                                  &
                           wgt*odifrad(2,klev,jsta) *                   &
                                  ABS(vazmrad(jsta))
                    END IF
                    IF(qualrad(3,klev,jsta) > 0) THEN
                      knt(5)=knt(5)+1
                      wgt=texp(-dist*rngsqi(5) -dzsq*zsqinv)
                      wgtsum(5)=wgtsum(5)+wgt
                      zsum(5)=zsum(5)+                                  &
                           wgt*odifrad(3,klev,jsta)
                    END IF
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval obs
!
!-----------------------------------------------------------------------
!
        IF(retsw > 0) THEN
          DO jsta=1,ncolret
            IF(iuseret(isrcret(iret(jsta))) > 0) THEN
              dist=(xua(ista)-xretc(jsta))*(xua(ista)-xretc(jsta))      &
                  +(yua(ista)-yretc(jsta))*(yua(ista)-yretc(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevret(jsta)
                  dzsq=(hgtua(k,ista)-hgtretc(klev,jsta)) *             &
                       (hgtua(k,ista)-hgtretc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualret(ivar,klev,jsta) > 0) THEN
                        knt(ivar)=knt(ivar)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar)=wgtsum(ivar) + wgt
                        zsum(ivar)=zsum(ivar)+                          &
                            wgt*odifret(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Divide by sum of weights
!
!-----------------------------------------------------------------------
!
        DO ivar=1,nvar
          IF(wgtsum(ivar) > wlim) oanxua(ivar,k,ista)=oanxua(ivar,k,ista) + &
              zsum(ivar)/wgtsum(ivar)
        END DO
      END DO
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop over all radar columns
!
!-----------------------------------------------------------------------
!
  DO ista=1,ncolrad
    IF(iuserad(isrcrad(irad(ista))) > 0) THEN
      DO k=1,nlevrad(ista)
        IF(qualrad(1,k,ista) > 0 .AND. qualrad(2,k,ista) > 0 ) THEN
          DO ivar=1,nvar
            zsum(ivar)=0.
            wgtsum(ivar)=0.
            knt(ivar)=0
          END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
          IF(sngsw > 0) THEN
            DO jsta=1,nobsng
             IF(isrcsng(jsta) > 0) THEN
              IF(iusesng(isrcsng(jsta)) > 0) THEN
                dist=(xradc(ista)-xsng(jsta))*(xradc(ista)-xsng(jsta))  &
                    +(yradc(ista)-ysng(jsta))*(yradc(ista)-ysng(jsta))
                IF(dist < rlimsq) THEN
                  dzsq=(hgtradc(k,ista)-hgtsng(jsta)) *                 &
                       (hgtradc(k,ista)-hgtsng(jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualsng(ivar,jsta) > 0) THEN
                        knt(ivar)=knt(ivar)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar)=wgtsum(ivar) + wgt
                        zsum(ivar)=zsum(ivar)+                          &
                             wgt*odifsng(ivar,jsta)
                      END IF
                    END DO
                  END IF
                END IF     ! horizontal dist check
              END IF     ! source check
             END IF      ! isrcsng check
            END DO
          END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
          IF(uasw > 0) THEN
            DO jsta=1,nobsua
              IF(iuseua(isrcua(jsta)) > 0) THEN
                dist=(xradc(ista)-xua(jsta))*(xradc(ista)-xua(jsta))    &
                    +(yradc(ista)-yua(jsta))*(yradc(ista)-yua(jsta))
                IF(dist < rlimsq) THEN
                  DO klev=1,nlevsua(jsta)
                    dzsq=(hgtradc(k,ista)-hgtua(klev,jsta)) *           &
                         (hgtradc(k,ista)-hgtua(klev,jsta))
                    IF(dzsq < zlimsq) THEN
                      DO ivar=1,nvar
                        IF(qualua(ivar,klev,jsta) > 0) THEN
                          knt(ivar)=knt(ivar)+1
                          wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                          wgtsum(ivar)=wgtsum(ivar) + wgt
                          zsum(ivar)=zsum(ivar)+                        &
                              wgt*odifua(ivar,klev,jsta)
                        END IF
                      END DO
                    END IF
                  END DO
                END IF     ! horizontal dist check
              END IF     ! source check
            END DO
          END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
          IF(radsw > 0) THEN
            DO jsta=1,ncolrad
              IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
                dist=(xradc(ista)-xradc(jsta))*(xradc(ista)-xradc(jsta)) &
                    +(yradc(ista)-yradc(jsta))*(yradc(ista)-yradc(jsta))
                IF(dist < rlimsq) THEN
                  DO klev=1,nlevrad(jsta)
                    dzsq=(hgtradc(k,ista)-hgtradc(klev,jsta)) *         &
                         (hgtradc(k,ista)-hgtradc(klev,jsta))
                    IF(dzsq < zlimsq) THEN
                      IF(qualrad(1,klev,jsta) > 0) THEN
                        knt(1)=knt(1)+1
                        wgt=texp(-dist*rngsqi(1) -dzsq*zsqinv)
                        wgtsum(1)=wgtsum(1)+wgt*                        &
                               ABS(uazmrad(jsta))
                        zsum(1)=zsum(1)+                                &
                            wgt*odifrad(1,klev,jsta) *                  &
                                   ABS(uazmrad(jsta))
                      END IF
                      IF(qualrad(2,klev,jsta) > 0) THEN
                        knt(2)=knt(2)+1
                        wgt=texp(-dist*rngsqi(2) -dzsq*zsqinv)
                        wgtsum(2)=wgtsum(2)+wgt*                        &
                               ABS(vazmrad(jsta))
                        zsum(2)=zsum(2)+                                &
                            wgt*odifrad(2,klev,jsta) *                  &
                                  ABS(vazmrad(jsta))
                      END IF
                      IF(qualrad(3,klev,jsta) > 0) THEN
                        knt(5)=knt(5)+1
                        wgt=texp(-dist*rngsqi(5) -dzsq*zsqinv)
                        wgtsum(5)=wgtsum(5)+wgt
                        zsum(5)=zsum(5)+                                &
                            wgt*odifrad(3,klev,jsta)
                      END IF
                    END IF
                  END DO
                END IF     ! horizontal dist check
              END IF     ! source check
            END DO
          END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar retrieval obs
!
!-----------------------------------------------------------------------
!
          IF(retsw > 0) THEN
            DO jsta=1,ncolret
              IF(iuseret(isrcret(iret(jsta))) > 0) THEN
                dist=(xradc(ista)-xretc(jsta))*(xradc(ista)-xretc(jsta)) &
                    +(yradc(ista)-yretc(jsta))*(yradc(ista)-yretc(jsta))
                IF(dist < rlimsq) THEN
                  DO klev=1,nlevret(jsta)
                    dzsq=(hgtradc(k,ista)-hgtretc(klev,jsta)) *         &
                         (hgtradc(k,ista)-hgtretc(klev,jsta))
                    IF(dzsq < zlimsq) THEN
                      DO ivar=1,nvar
                        IF(qualret(ivar,klev,jsta) > 0) THEN
                          knt(ivar)=knt(ivar)+1
                          wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                          wgtsum(ivar)=wgtsum(ivar) + wgt
                          zsum(ivar)=zsum(ivar)+                        &
                               wgt*odifret(ivar,klev,jsta)
                        END IF
                      END DO
                    END IF
                  END DO
                END IF     ! horizontal dist check
              END IF     ! source check
            END DO
          END IF
!
!-----------------------------------------------------------------------
!
!  Divide by sum of weights
!
!-----------------------------------------------------------------------
!
          IF(wgtsum(1) > wlim) oanxrad(1,k,ista)=oanxrad(1,k,ista) +    &
                zsum(1)/wgtsum(1)
          IF(wgtsum(2) > wlim) oanxrad(2,k,ista)=oanxrad(2,k,ista) +    &
                zsum(2)/wgtsum(2)
          IF(wgtsum(5) > wlim) oanxrad(3,k,ista)=oanxrad(3,k,ista) +    &
                zsum(5)/wgtsum(5)
        END IF
      END DO
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop over all radar retrievals
!
!-----------------------------------------------------------------------
!
  DO ista=1,ncolret
    IF(iuseret(isrcret(iret(ista))) > 0) THEN
      DO k=1,nlevret(ista)
        DO ivar=1,nvar
          zsum(ivar)=0.
          wgtsum(ivar)=0.
          knt(ivar)=0
        END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
        IF(sngsw > 0) THEN
          DO jsta=1,nobsng
           IF(isrcsng(jsta) > 0) THEN
            IF(iusesng(isrcsng(jsta)) > 0) THEN
              dist=(xretc(ista)-xsng(jsta))*(xretc(ista)-xsng(jsta))    &
                  +(yretc(ista)-ysng(jsta))*(yretc(ista)-ysng(jsta))
              IF(dist < rlimsq) THEN
                dzsq=(hgtretc(k,ista)-hgtsng(jsta)) *                   &
                     (hgtretc(k,ista)-hgtsng(jsta))
                IF(dzsq < zlimsq) THEN
                  DO ivar=1,nvar
                    IF(qualsng(ivar,jsta) > 0) THEN
                      knt(ivar)=knt(ivar)+1
                      wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                      wgtsum(ivar)=wgtsum(ivar) + wgt
                      zsum(ivar)=zsum(ivar)+                            &
                           wgt*odifsng(ivar,jsta)
                    END IF
                  END DO
                END IF
              END IF     ! horizontal dist check
            END IF     ! source check
           END IF      ! isrcsng check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
        IF(uasw > 0) THEN
          DO jsta=1,nobsua
            IF(iuseua(isrcua(jsta)) > 0) THEN
              dist=(xretc(ista)-xua(jsta))*(xretc(ista)-xua(jsta))      &
                  +(yretc(ista)-yua(jsta))*(yretc(ista)-yua(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevsua(jsta)
                  dzsq=(hgtretc(k,ista)-hgtua(klev,jsta)) *             &
                       (hgtretc(k,ista)-hgtua(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualua(ivar,klev,jsta) > 0) THEN
                        knt(ivar)=knt(ivar)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar)=wgtsum(ivar) + wgt
                        zsum(ivar)=zsum(ivar)+                          &
                            wgt*odifua(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
        IF(radsw > 0) THEN
          DO jsta=1,ncolrad
            IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
              dist=(xretc(ista)-xradc(jsta))*(xretc(ista)-xradc(jsta))  &
                  +(yretc(ista)-yradc(jsta))*(yretc(ista)-yradc(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevrad(jsta)
                  dzsq=(hgtretc(k,ista)-hgtradc(klev,jsta)) *           &
                       (hgtretc(k,ista)-hgtradc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    IF(qualrad(1,klev,jsta) > 0) THEN
                      knt(1)=knt(1)+1
                      wgt=texp(-dist*rngsqi(1) -dzsq*zsqinv)
                      wgtsum(1)=wgtsum(1)+wgt*                          &
                               ABS(uazmrad(jsta))
                      zsum(1)=zsum(1)+                                  &
                           wgt*odifrad(1,klev,jsta) *                   &
                               ABS(uazmrad(jsta))
                    END IF
                    IF(qualrad(2,klev,jsta) > 0) THEN
                      knt(2)=knt(2)+1
                      wgt=texp(-dist*rngsqi(2) -dzsq*zsqinv)
                      wgtsum(2)=wgtsum(2)+wgt*                          &
                               ABS(vazmrad(jsta))
                      zsum(2)=zsum(2)+                                  &
                           wgt*odifrad(2,klev,jsta) *                   &
                                  ABS(vazmrad(jsta))
                    END IF
                    IF(qualrad(3,klev,jsta) > 0) THEN
                      knt(5)=knt(5)+1
                      wgt=texp(-dist*rngsqi(5) -dzsq*zsqinv)
                      wgtsum(5)=wgtsum(5)+wgt
                      zsum(5)=zsum(5)+                                  &
                           wgt*odifrad(3,klev,jsta)
                    END IF
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval obs
!
!-----------------------------------------------------------------------
!
        IF(retsw > 0) THEN
          DO jsta=1,ncolret
            IF(iuseret(isrcret(jsta)) > 0) THEN
              dist=(xretc(ista)-xretc(jsta))*(xretc(ista)-xretc(jsta))  &
                  +(yretc(ista)-yretc(jsta))*(yretc(ista)-yretc(jsta))
              IF(dist < rlimsq) THEN
                DO klev=1,nlevret(jsta)
                  dzsq=(hgtretc(k,ista)-hgtretc(klev,jsta)) *           &
                       (hgtretc(k,ista)-hgtretc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualret(ivar,klev,jsta) > 0) THEN
                        knt(ivar)=knt(ivar)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar)=wgtsum(ivar) + wgt
                        zsum(ivar)=zsum(ivar)+                          &
                            wgt*odifret(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END IF     ! horizontal dist check
            END IF     ! source check
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
!  Divide by sum of weights
!
!-----------------------------------------------------------------------
!
        DO ivar=1,nvar
          IF(wgtsum(ivar) > wlim) oanxret(ivar,k,ista)=oanxret(ivar,k,ista) + &
              zsum(ivar)/wgtsum(ivar)
        END DO
      END DO
    END IF
  END DO
!
  istatus=0
  RETURN
END SUBROUTINE barobs3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BARGRD3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bargrd3d(nx,ny,nz,                                           &
           nvar,nvarrad,nzua,nzrdr,nzret,                               &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,xs,ys,zs,                     &
           odifsng,xsng,ysng,hgtsng,qualsng,isrcsng,nobsng,             &
           odifua,xua,yua,hgtua,qualua,isrcua,nlevsua,nobsua,           &
           odifrad,xradc,yradc,hgtradc,                                 &
           uazmrad,vazmrad,qualrad,irad,isrcrad,nlevrad,ncolrad,        &
           odifret,xretc,yretc,hgtretc,                                 &
           qualret,iret,isrcret,nlevret,ncolret,                        &
           iusesng,iuseua,iuserad,iuseret,                              &
           sngsw,uasw,radsw,retsw,kpvrsq,rpass,                         &
           wlim,zrngsq,zwlim,                                           &
           rngsqi,knt,wgtsum,zsum,                                      &
           anx,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Objectively analyze data on a regular grid specified
!  by x and y.  Barnes weight function is employed.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  September, 1995
!
!  MODIFICATION HISTORY:
!
!  July, 1996 (KB)
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!
!  Sizing Arguments
!
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nvarrad,nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret
!
!-----------------------------------------------------------------------
!
!  Grid arguments
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: odifsng(nvar,mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: nobsng
  REAL :: odifua(nvar,nzua,mxua)             ! variable to analyse
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Scratch space
!
!-----------------------------------------------------------------------
!
  REAL :: rngsqi(nvar)
  INTEGER :: knt(nvar,nz)
  REAL :: wgtsum(nvar,nz)
  REAL :: zsum(nvar,nz)
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)
  INTEGER :: sngsw,uasw,radsw,retsw
  REAL :: kpvrsq(nvar)
  REAL :: rpass,                                                        &
       wlim
  REAL :: zrngsq,                                                       &
       zwlim
INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Output grid arguments
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL :: texp
  REAL :: rlimsq,zlimsq,zsqinv
  REAL :: dist,dzsq,wgt
  INTEGER :: i,j,k,jsta,klev,ivar
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rlimsq=-rpass*ALOG(wlim)
  DO ivar=1,nvar
    rngsqi(ivar)=1./(kpvrsq(ivar)*rpass)
  END DO
!
  zlimsq=-zrngsq*ALOG(zwlim)
  zsqinv=1./zrngsq
!
!-----------------------------------------------------------------------
!
!  Loop over all grid points
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      DO k=1,nz-1
        DO ivar=1,nvar
          knt(ivar,k)=0
          wgtsum(ivar,k)=0.
          zsum(ivar,k)=0.
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
      IF(sngsw > 0) THEN
        DO jsta=1,nobsng
         IF(isrcsng(jsta) > 0) THEN
          IF(iusesng(isrcsng(jsta)) > 0) THEN
            dist=(xs(i)-xsng(jsta))*(xs(i)-xsng(jsta))                  &
                +(ys(j)-ysng(jsta))*(ys(j)-ysng(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nz-1
                dzsq=(zs(i,j,k)-hgtsng(jsta)) *                         &
                     (zs(i,j,k)-hgtsng(jsta))
                IF(dzsq < zlimsq) THEN
                  DO ivar=1,nvar
                    IF(qualsng(ivar,jsta) > 0) THEN
                      knt(ivar,k)=knt(ivar,k)+1
                      wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                      wgtsum(ivar,k)=wgtsum(ivar,k) + wgt
                      zsum(ivar,k)=zsum(ivar,k)+                        &
                            wgt*odifsng(ivar,jsta)
                    END IF
                  END DO
                END IF
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
         END IF      ! isrcsng check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
      IF(uasw > 0) THEN
        DO jsta=1,nobsua
          IF(iuseua(isrcua(jsta)) > 0) THEN
            dist=(xs(i)-xua(jsta))*(xs(i)-xua(jsta))                    &
                +(ys(j)-yua(jsta))*(ys(j)-yua(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nz-1
                DO klev=1,nlevsua(jsta)
                  dzsq=(zs(i,j,k)-hgtua(klev,jsta)) *                   &
                       (zs(i,j,k)-hgtua(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualua(ivar,klev,jsta) > 0) THEN
                        knt(ivar,k)=knt(ivar,k)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar,k)=wgtsum(ivar,k) + wgt
                        zsum(ivar,k)=zsum(ivar,k)+                      &
                               wgt*odifua(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
      IF(radsw > 0) THEN
        DO jsta=1,ncolrad
          IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
            dist=(xs(i)-xradc(jsta))*(xs(i)-xradc(jsta))                &
                +(ys(j)-yradc(jsta))*(ys(j)-yradc(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nz-1
                DO klev=1,nlevrad(jsta)
                  dzsq=(zs(i,j,k)-hgtradc(klev,jsta)) *                 &
                       (zs(i,j,k)-hgtradc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    IF(qualrad(1,klev,jsta) > 0) THEN
                      knt(1,k)=knt(1,k)+1
                      wgt=texp(-dist*rngsqi(1) -dzsq*zsqinv)
                      wgtsum(1,k)=wgtsum(1,k)+                          &
                                  wgt*uazmrad(jsta)
                      zsum(1,k)=zsum(1,k)+                              &
                          wgt*odifrad(1,klev,jsta) *                    &
                                   uazmrad(jsta)
                    END IF
                    IF(qualrad(2,klev,jsta) > 0) THEN
                      knt(2,k)=knt(2,k)+1
                      wgt=texp(-dist*rngsqi(2) -dzsq*zsqinv)
                      wgtsum(2,k)=wgtsum(2,k)+                          &
                                  wgt*vazmrad(jsta)
                      zsum(2,k)=zsum(2,k)+                              &
                          wgt*odifrad(2,klev,jsta) *                    &
                                   vazmrad(jsta)
                    END IF
                    IF(qualrad(3,klev,jsta) > 0) THEN
                      knt(5,k)=knt(5,k)+1
                      wgt=texp(-dist*rngsqi(5) -dzsq*zsqinv)
                      wgtsum(5,k)=wgtsum(5,k)+wgt
                      zsum(5,k)=zsum(5,k)+                              &
                          wgt*odifrad(3,klev,jsta)
                    END IF
                  END IF
                END DO
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval obs
!
!-----------------------------------------------------------------------
!
      IF(retsw > 0) THEN
        DO jsta=1,ncolret
          IF(iuseret(isrcret(iret(jsta))) > 0) THEN
            dist=(xs(i)-xretc(jsta))*(xs(i)-xretc(jsta))                &
                +(ys(j)-yretc(jsta))*(ys(j)-yretc(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nz-1
                DO klev=1,nlevret(jsta)
                  dzsq=(zs(i,j,k)-hgtretc(klev,jsta)) *                 &
                       (zs(i,j,k)-hgtretc(klev,jsta))
                  IF(dzsq < zlimsq) THEN
                    DO ivar=1,nvar
                      IF(qualret(ivar,klev,jsta) > 0) THEN
                        knt(ivar,k)=knt(ivar,k)+1
                        wgt=texp(-dist*rngsqi(ivar) -dzsq*zsqinv)
                        wgtsum(ivar,k)=wgtsum(ivar,k) + wgt
                        zsum(ivar,k)=zsum(ivar,k)+                      &
                               wgt*odifret(ivar,klev,jsta)
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END IF     ! horizontal dist check
          END IF     ! source check
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Divide by sum of weights
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO ivar=1,nvar
          IF(wgtsum(ivar,k) > wlim) anx(i,j,k,ivar)=anx(i,j,k,ivar) +   &
                            zsum(ivar,k)/wgtsum(ivar,k)
        END DO
      END DO
    END DO
  END DO
!
  istatus=0
  RETURN
END SUBROUTINE bargrd3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE OBSCOR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE obscor(nvar,nvarrad,nzua,nzrdr,nzret,                        &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                    &
           qobsng,xsng,ysng,hgtsng,trnsng,                              &
           qualsng,isrcsng,icatsng,nobsng,                              &
           qobsua,xua,yua,hgtua,trnua,                                  &
           qualua,isrcua,nlevsua,nobsua,                                &
           qobsrad,xradc,yradc,hgtradc,trnradc,                         &
           uazmrad,vazmrad,                                             &
           qualrad,irad,isrcrad,nlevrad,ncolrad,                        &
           qobsret,xretc,yretc,hgtretc,trnretc,                         &
           qualret,iret,isrcret,nlevret,ncolret,                        &
           iusesng,iuseua,iuserad,iuseret,                              &
           sngsw,uasw,radsw,retsw,kpvrsq,rpass,rngsqi,                  &
           wlim,zrngsq,trnropt,trnrcst,trnrngsq,                        &
           corsng,corua,corrad,corret,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate covariance sums at observation points for
!  Bratseth analysis method.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (KB)
!  Added radar data and other improvements.
!
!  July, 1996 (KB)
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!
!  Jan 5, 1997 (KB)
!  New version to reorder summing to speed-up code.
!
!  Oct 30, 2006 (KWT)
!  Code fails to check if some observations have been marked to not be
!  used.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Sizing arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nvar,nvarrad,nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  REAL :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: qobsng(nvar,mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: trnsng(mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
!
  REAL :: qobsua(nvar,nzua,mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: trnua(mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL :: qobsrad(nvarrad,nzrdr,mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: trnradc(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
  REAL :: qobsret(nvar,nzret,mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: trnretc(mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)
  INTEGER :: sngsw,uasw,radsw,retsw
  REAL :: kpvrsq(nvar)
  REAL :: rpass,                                                        &
       wlim
  REAL :: zrngsq
  REAL :: rngsqi(nvar)
  INTEGER :: trnropt
  REAL :: trnrcst
  REAL :: trnrngsq
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  REAL :: corsng(nvar,mxsng)
  REAL :: corua(nvar,nzua,mxua)
  REAL :: corrad(nvarrad,nzrdr,mxcolrad)
  REAL :: corret(nvar,nzret,mxcolret)
  INTEGER :: istatus
!
  REAL :: ftabinv
  INTEGER :: ntabexp
  PARAMETER (ntabexp=5000)
  REAL :: tabexp(ntabexp)
  COMMON /ftabexp/ ftabinv
  COMMON /tablexp/ tabexp
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL :: rsq,rlim,rlimsq,zsqinv,trnsqinv
  REAL :: dist,dzsq,dtrn2,dtrnsq,rnorm
  REAL :: aglista,agljsta,aglavg
  INTEGER :: ista,jsta,k,klev,ivar,indextab

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rlimsq=0.
  DO ivar=1,nvar
    rsq=(kpvrsq(ivar)*rpass)
    rngsqi(ivar)=1./rsq
    rlimsq=AMAX1(rlimsq,rsq)
  END DO
  rlimsq=-rlimsq*ALOG(wlim)
  rlim=0.001*SQRT(rlimsq)
  IF (myproc == 0) WRITE(6,'(/a,f10.2,a)')                              &
       ' OBSCOR: Influence cutoff radius: ',rlim,' km.'

  zsqinv=1./zrngsq
!
  IF(trnropt > 0) THEN
    trnsqinv=1./trnrngsq
  ELSE
    trnsqinv=0.
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize covariance sums with observation error variance
!
!-----------------------------------------------------------------------
!
! OpenMP:
!$OMP PARALLEL DO PRIVATE(ista,ivar)
  DO ista=1,nobsng
    DO ivar=1,nvar
      corsng(ivar,ista)=qobsng(ivar,ista)
    END DO
  END DO
!$OMP PARALLEL DO PRIVATE(ista,k,ivar)
  DO ista=1,nobsua
    DO k=1,nlevsua(ista)
      DO ivar=1,nvar
        corua(ivar,k,ista)=qobsua(ivar,k,ista)
      END DO
    END DO
  END DO
!$OMP PARALLEL DO PRIVATE(ista,k,ivar)
  DO ista=1,ncolret
    DO k=1,nlevret(ista)
      DO ivar=1,nvar
        corret(ivar,k,ista)=qobsret(ivar,k,ista)
      END DO
    END DO
  END DO
!$OMP PARALLEL DO PRIVATE(ista,k,ivar)
  DO ista=1,ncolrad
    DO k=1,nlevrad(ista)
      DO ivar=1,nvarrad
        corrad(ivar,k,ista)=qobsrad(ivar,k,ista)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
  IF(sngsw > 0) THEN
    DO jsta=1,nobsng
     IF(isrcsng(jsta) > 0) THEN
      IF(iusesng(isrcsng(jsta)) > 0) THEN
        agljsta=max(0.,(hgtsng(jsta)-trnsng(jsta)))
!
!-----------------------------------------------------------------------
!
!  Add to single-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xsng(jsta))*(xsng(ista)-xsng(jsta))          &
                +(ysng(ista)-ysng(jsta))*(ysng(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dzsq=zsqinv*(hgtsng(ista)-hgtsng(jsta)) *                   &
                          (hgtsng(ista)-hgtsng(jsta))
              aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
              aglavg=(0.5*(agljsta+aglista))+trnrcst
              dtrnsq=trnsqinv*(trnsng(ista)-trnsng(jsta))*                &
                         (trnsng(ista)-trnsng(jsta))/(aglavg*aglavg)
              DO ivar=1,nvar
                IF(qualsng(ivar,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corsng(ivar,ista)=corsng(ivar,ista) +                   &
                      tabexp(indextab)*xcor(icatsng(ista),icatsng(jsta))
                END IF
              END DO
            END IF
          END IF
         END IF      ! isrcsng check
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xsng(jsta))*(xua(ista)-xsng(jsta))            &
                +(yua(ista)-ysng(jsta))*(yua(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnua(ista)-trnsng(jsta))*                           &
                    (trnua(ista)-trnsng(jsta))
              DO k=1,nlevsua(ista)
                dzsq=zsqinv*(hgtua(k,ista)-hgtsng(jsta)) *                &
                             (hgtua(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                DO ivar=1,nvar
                  IF(qualsng(ivar,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corua(ivar,k,ista)=corua(ivar,k,ista)+tabexp(indextab)
                  END IF
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to retrieval covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xsng(jsta))*(xretc(ista)-xsng(jsta))        &
                +(yretc(ista)-ysng(jsta))*(yretc(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnsng(jsta))*                         &
                    (trnretc(ista)-trnsng(jsta))
              DO k=1,nlevret(ista)
                dzsq=zsqinv*(hgtretc(k,ista)-hgtsng(jsta)) *              &
                            (hgtretc(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                DO ivar=1,nvar
                  IF(qualsng(ivar,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corret(ivar,k,ista)=                                  &
                        corret(ivar,k,ista)+tabexp(indextab)
                  END IF
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar data covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xsng(jsta))*(xradc(ista)-xsng(jsta))        &
                +(yradc(ista)-ysng(jsta))*(yradc(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnsng(jsta))*                         &
                    (trnradc(ista)-trnsng(jsta))
              DO k=1,nlevrad(ista)
                dzsq=zsqinv*(hgtradc(k,ista)-hgtsng(jsta)) *              &
                     (hgtradc(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                IF(qualsng(1,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corrad(1,k,ista)=                                       &
                      corrad(1,k,ista)+tabexp(indextab)
                END IF
!
                IF(qualsng(2,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corrad(2,k,ista)=                                       &
                      corrad(2,k,ista)+tabexp(indextab)
                END IF
!
                IF(qualsng(5,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corrad(3,k,ista)=                                       &
                      corrad(3,k,ista)+tabexp(indextab)
                END IF
              END DO
            END IF
          END IF
        END DO
      END IF     ! source check
     END IF      ! isrcsng check
    END DO
  END IF   ! sngsw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
  IF(uasw > 0) THEN
    DO jsta=1,nobsua
      IF(iuseua(isrcua(jsta)) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xua(jsta))*(xsng(ista)-xua(jsta))            &
                +(ysng(ista)-yua(jsta))*(ysng(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtua(klev,jsta)) *             &
                            (hgtsng(ista)-hgtua(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnua(jsta))*               &
                           (trnsng(ista)-trnua(jsta))/(aglavg*aglavg)
                DO ivar=1,nvar
                  IF(qualua(ivar,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corsng(ivar,ista)=                                    &
                        corsng(ivar,ista)+tabexp(indextab)
                  END IF
                END DO ! ivar
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
           IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xua(jsta))*(xua(ista)-xua(jsta))              &
                +(yua(ista)-yua(jsta))*(yua(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnua(ista)-trnua(jsta))*                            &
                    (trnua(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtua(klev,jsta)) *          &
                              (hgtua(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                  DO ivar=1,nvar
                    IF(qualua(ivar,klev,jsta) > 0) THEN
                      rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                      indextab=min((nint(rnorm)+1),ntabexp)
                      corua(ivar,k,ista)=                                 &
                          corua(ivar,k,ista)+tabexp(indextab)
                    END IF
                  END DO  ! ivar
                END DO  ! k
              END DO  ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to retrieval covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xua(jsta))*(xretc(ista)-xua(jsta))          &
                +(yretc(ista)-yua(jsta))*(yretc(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnua(jsta))*                          &
                    (trnretc(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtua(klev,jsta)) *        &
                              (hgtretc(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                  DO ivar=1,nvar
                    IF(qualua(ivar,klev,jsta) > 0) THEN
                      rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                      indextab=min((nint(rnorm)+1),ntabexp)
                      corret(ivar,k,ista)=                                &
                          corret(ivar,k,ista)+tabexp(indextab)
                    END IF
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xua(jsta))*(xradc(ista)-xua(jsta))          &
                +(yradc(ista)-yua(jsta))*(yradc(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnua(jsta))*                          &
                    (trnradc(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtua(klev,jsta)) *        &
                              (hgtradc(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  IF(qualua(1,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(1,k,ista)=                                    &
                        corrad(1,k,ista)+tabexp(indextab)
                  END IF
!
                  IF(qualua(2,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(2,k,ista)=                                    &
                        corrad(2,k,ista)+tabexp(indextab)
                  END IF
!
                  IF(qualua(5,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(3,k,ista)=                                    &
                        corrad(3,k,ista)+tabexp(indextab)
                  END IF
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO
      END IF     ! source check
    END DO
  END IF   ! uasw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval data
!
!-----------------------------------------------------------------------
!
  IF(retsw > 0) THEN
    DO jsta=1,ncolret
      IF(iuseret(isrcret(iret(jsta))) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xretc(jsta))*(xsng(ista)-xretc(jsta))        &
                +(ysng(ista)-yretc(jsta))*(ysng(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtretc(klev,jsta)) *           &
                            (hgtsng(ista)-hgtretc(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnretc(jsta))*             &
                           (trnsng(ista)-trnretc(jsta))/(aglavg*aglavg)
                DO ivar=1,nvar
                  IF(qualret(ivar,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corsng(ivar,ista)=                &
                        corsng(ivar,ista)+tabexp(indextab)
                  END IF
                END DO ! ivar
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xretc(jsta))*(xua(ista)-xretc(jsta))          &
                +(yua(ista)-yretc(jsta))*(yua(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnua(jsta))*                          &
                    (trnretc(ista)-trnua(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtretc(klev,jsta)) *        &
                              (hgtua(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                  DO ivar=1,nvar
                    IF(qualret(ivar,klev,jsta) > 0 ) THEN
                      rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                      indextab=min((nint(rnorm)+1),ntabexp)
                      corua(ivar,k,ista)=                                &
                          corua(ivar,k,ista)+tabexp(indextab)
                    END IF
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to retrieval covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xretc(jsta))*(xretc(ista)-xretc(jsta))      &
                +(yretc(ista)-yretc(jsta))*(yretc(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnretc(jsta))*                        &
                    (trnretc(ista)-trnretc(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtretc(klev,jsta)) *      &
                              (hgtretc(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
                  DO ivar=1,nvar
                    IF(qualret(ivar,klev,jsta) > 0 ) THEN
                      rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                      indextab=min((nint(rnorm)+1),ntabexp)
                      corret(ivar,k,ista)=                               &
                          corret(ivar,k,ista)+tabexp(indextab)
                    END IF
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to radar covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xretc(jsta))*(xradc(ista)-xretc(jsta))      &
                +(yradc(ista)-yretc(jsta))*(yradc(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnretc(jsta))*                        &
                    (trnradc(ista)-trnretc(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtretc(klev,jsta)) *      &
                              (hgtradc(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  IF(qualret(1,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(1,k,ista)=                                     &
                        corrad(1,k,ista)+tabexp(indextab)
                  END IF
!
                  IF(qualret(2,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(2,k,ista)=                                     &
                        corrad(2,k,ista)+tabexp(indextab)
                  END IF
!
                  IF(qualret(5,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(3,k,ista)=                                     &
                        corrad(3,k,ista)+tabexp(indextab)
                  END IF
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
      END IF ! source check
    END DO
  END IF ! retsw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
  IF(radsw > 0) THEN
    DO jsta=1,ncolrad
      IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xradc(jsta))*(xsng(ista)-xradc(jsta))        &
                +(ysng(ista)-yradc(jsta))*(ysng(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtradc(klev,jsta)) *           &
                            (hgtsng(ista)-hgtradc(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnradc(jsta))*             &
                       (trnsng(ista)-trnradc(jsta))/(aglavg*aglavg)
!
                IF(qualrad(1,klev,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corsng(1,ista)=corsng(1,ista)+                          &
                               (tabexp(indextab)*ABS(uazmrad(jsta)))
                END IF
!
                IF(qualrad(2,klev,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corsng(2,ista)=corsng(2,ista)+                          &
                               (tabexp(indextab)*ABS(vazmrad(jsta)))
                END IF
!
                IF(qualrad(3,klev,jsta) > 0) THEN
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  corsng(5,ista)=corsng(5,ista)+                          &
                               tabexp(indextab)
                END IF
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xradc(jsta))*(xua(ista)-xradc(jsta))          &
                +(yua(ista)-yradc(jsta))*(yua(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnua(ista)-trnradc(jsta))*                          &
                    (trnua(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtradc(klev,jsta)) *        &
                              (hgtua(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  IF(qualrad(1,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corua(1,k,ista)=corua(1,k,ista)+                      &
                               (tabexp(indextab)*ABS(uazmrad(jsta)))
                  END IF
!
                  IF(qualrad(2,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corua(2,k,ista)=corua(2,k,ista)+                      &
                               (tabexp(indextab)*ABS(vazmrad(jsta)))
                  END IF
!
                  IF(qualrad(3,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corua(5,k,ista)=corua(5,k,ista)+tabexp(indextab)
                  END IF
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to retrieval covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xradc(jsta))*(xretc(ista)-xradc(jsta))      &
                +(yretc(ista)-yradc(jsta))*(yretc(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnradc(jsta))*                        &
                    (trnretc(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtradc(klev,jsta)) *      &
                              (hgtretc(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  IF(qualrad(1,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corret(1,k,ista)=corret(1,k,ista)+&
                               (tabexp(indextab)*ABS(uazmrad(jsta)))
                  END IF
!
                  IF(qualrad(2,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corret(2,k,ista)=corret(2,k,ista) +                  &
                               (tabexp(indextab)*ABS(vazmrad(jsta)))
                  END IF
!
                  IF(qualrad(3,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corret(5,k,ista)=                                    &
                        corret(5,k,ista)+tabexp(indextab)
                  END IF
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar covariance sum
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xradc(jsta))*(xradc(ista)-xradc(jsta))      &
                +(yradc(ista)-yradc(jsta))*(yradc(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnradc(jsta))*                        &
                    (trnradc(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtradc(klev,jsta)) *      &
                              (hgtradc(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  IF(qualrad(1,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(1,k,ista)=corrad(1,k,ista) +                   &
                               (tabexp(indextab)*ABS(uazmrad(jsta)))
                  END IF
!
                  IF(qualrad(2,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(2,k,ista)=corrad(2,k,ista) +                   &
                               (tabexp(indextab)*ABS(vazmrad(jsta)))
                  END IF
!
                  IF(qualrad(3,klev,jsta) > 0) THEN
                    rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    corrad(3,k,ista)=                                     &
                        corrad(3,k,ista)+tabexp(indextab)
                  END IF
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
      END IF     ! source check
    END DO
  END IF   ! radsw check
!
!-----------------------------------------------------------------------
!
!  Take inverse to speed up later calculations or
!  set covariance inverse to zero as a mask for "bad" obs
!
!-----------------------------------------------------------------------
!
  DO ista=1,nobsng
    DO ivar=1,nvar
      IF(qualsng(ivar,ista) > 0 ) THEN
        corsng(ivar,ista)=1./corsng(ivar,ista)
      ELSE
        corsng(ivar,ista)=0.
      END IF
    END DO
  END DO
!
  DO ista=1,nobsua
    DO k=1,nlevsua(ista)
      DO ivar=1,nvar
        IF(qualua(ivar,k,ista) > 0 ) THEN
          corua(ivar,k,ista)=1./corua(ivar,k,ista)
        ELSE
          corua(ivar,k,ista)=0.
        END IF
      END DO
    END DO
  END DO
!
  DO ista=1,ncolret
    DO k=1,nlevret(ista)
      DO ivar=1,nvar
        IF(qualret(ivar,k,ista) > 0 ) THEN
          corret(ivar,k,ista)=1./corret(ivar,k,ista)
        ELSE
          corret(ivar,k,ista)=0.
        END IF
      END DO
    END DO
  END DO
!
  DO ista=1,ncolrad
    DO k=1,nlevrad(ista)
      DO ivar=1,nvarrad
        IF(qualrad(ivar,k,ista) > 0 ) THEN
          IF(corrad(ivar,k,ista) == 0.) THEN
            PRINT *, ' zero corrad: ',ivar,k,ista
            PRINT *, ' qobsrad: ',qobsrad(ivar,k,ista)
            CALL arpsstop("obscor: zero corrad",1)
          END IF
          corrad(ivar,k,ista)=1./corrad(ivar,k,ista)
        ELSE
          corrad(ivar,k,ista)=0.
        END IF
      END DO
    END DO
  END DO
!
  istatus=0
  RETURN
END SUBROUTINE obscor
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRTOBS3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE brtobs3d(nvar,nvarrad,nzua,nzrdr,nzret,                      &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xcor,                    &
           odifsng,qobsng,corsng,xsng,ysng,hgtsng,trnsng,               &
           isrcsng,icatsng,nobsng,                                      &
           odifua,qobsua,corua,xua,yua,hgtua,trnua,                     &
           isrcua,nlevsua,nobsua,                                       &
           odifrad,qobsrad,corrad,xradc,yradc,hgtradc,trnradc,          &
           uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,                &
           odifret,qobsret,corret,xretc,yretc,hgtretc,trnretc,          &
           iret,isrcret,nlevret,ncolret,                                &
           iusesng,iuseua,iuserad,iuseret,                              &
           sngsw,uasw,radsw,retsw,kpvrsq,rpass,                         &
           wlim,zrngsq,trnropt,trnrcst,trnrngsq,rngsqi,                &
           oanxsng,oanxua,oanxrad,oanxret,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Objectively analyze data at observation pts
!  Bratseth scheme is employed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (KB)
!  Added radar data and other improvements.
!  Keith Brewster, CAPS
!
!  July, 1996 (KB)
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!
!  Jan 5, 1997 (KB)
!  New version to reorder summing to speed-up code.
!  Renamed to BRTGRD3D from BRTGRDZ3d and BRTGRDTH3d.
!  When called to analyze using delta-theta, potential
!  temperature data is passed into the height arrays.
!
!  Oct 30, 2006 (KWT)
!  Code fails to check if some observations have been marked to not be
!  used.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Sizing arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nvar,nvarrad,nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  REAL :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: odifsng(nvar,mxsng)
  REAL :: qobsng(nvar,mxsng)
  REAL :: corsng(nvar,mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: trnsng(mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
!
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: qobsua(nvar,nzua,mxua)
  REAL :: corua(nvar,nzua,mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: trnua(mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: qobsrad(nvarrad,nzrdr,mxcolrad)
  REAL :: corrad(nvarrad,nzrdr,mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: trnradc(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: qobsret(nvar,nzret,mxcolret)
  REAL :: corret(nvar,nzret,mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: trnretc(mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Scratch space
!
!-----------------------------------------------------------------------
!
  REAL :: rngsqi(nvar)
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)
  INTEGER :: sngsw,uasw,radsw,retsw
  REAL :: kpvrsq(nvar)
  REAL :: rpass
  REAL :: wlim
  REAL :: zrngsq
  INTEGER :: trnropt
  REAL :: trnrcst
  REAL :: trnrngsq
!
!-----------------------------------------------------------------------
!
!  Output arguments
!
!-----------------------------------------------------------------------
!
  REAL :: oanxsng(nvar,mxsng)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxret(nvar,nzret,mxcolret)
  INTEGER :: istatus
!
  REAL :: ftabinv
  INTEGER :: ntabexp
  PARAMETER (ntabexp=5000)
  REAL :: tabexp(ntabexp)
  COMMON /ftabexp/ ftabinv
  COMMON /tablexp/ tabexp

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: rsq,rlim,rlimsq,zsqinv,trnsqinv
  REAL :: dist,dzsq,dtrn2,dtrnsq,rnorm
  REAL :: aglista,agljsta,aglavg
  INTEGER :: ista,jsta,k,klev,ivar,indextab
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rlimsq=0.
  DO ivar=1,nvar
    rsq=(kpvrsq(ivar)*rpass)
    rngsqi(ivar)=1./rsq
    rlimsq=AMAX1(rlimsq,rsq)
  END DO
  rlimsq=-rlimsq*ALOG(wlim)
  rlim=0.001*SQRT(rlimsq)
  IF (myproc == 0)                                                       &
    WRITE(6,'(a,f10.2,a)') ' Influence cutoff radius: ',rlim,' km.'
!
  zsqinv=1./zrngsq
!
  IF(trnropt > 0) THEN
    trnsqinv=1./trnrngsq
  ELSE
    trnsqinv=0.
  END IF
!
!-----------------------------------------------------------------------
!
!  Sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
  IF(sngsw > 0) THEN
    DO jsta=1,nobsng
     IF(isrcsng(jsta) > 0) THEN
      IF(iusesng(isrcsng(jsta)) > 0) THEN
        agljsta=max(0.,(hgtsng(jsta)-trnsng(jsta)))
!
!-----------------------------------------------------------------------
!
!  Add to single-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xsng(jsta))*(xsng(ista)-xsng(jsta))          &
                +(ysng(ista)-ysng(jsta))*(ysng(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              dzsq=zsqinv*(hgtsng(ista)-hgtsng(jsta))*                    &
                          (hgtsng(ista)-hgtsng(jsta))
              aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
              aglavg=(0.5*(agljsta+aglista))+trnrcst
              dtrnsq=trnsqinv*(trnsng(ista)-trnsng(jsta))*       &
                     (trnsng(ista)-trnsng(jsta))/(aglavg*aglavg)

              DO ivar=1,nvar
                rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxsng(ivar,ista)=oanxsng(ivar,ista) +                   &
                    (tabexp(indextab)*xcor(icatsng(ista),icatsng(jsta))   &
                     *odifsng(ivar,jsta)*corsng(ivar,jsta))
              END DO
            END IF
          END IF
         END IF      ! isrcsng check
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xsng(jsta))*(xua(ista)-xsng(jsta))            &
                +(yua(ista)-ysng(jsta))*(yua(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nlevsua(ista)
                dzsq=zsqinv*(hgtua(k,ista)-hgtsng(jsta)) *                &
                            (hgtua(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnua(ista)-trnsng(jsta))*      &
                       (trnua(ista)-trnsng(jsta))/(aglavg*aglavg)
                DO ivar=1,nvar
                  rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxua(ivar,k,ista)=oanxua(ivar,k,ista) +               &
                   (tabexp(indextab)*odifsng(ivar,jsta)*corsng(ivar,jsta))
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to retrieval data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xsng(jsta))*(xretc(ista)-xsng(jsta))        &
                +(yretc(ista)-ysng(jsta))*(yretc(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nlevret(ista)
                dzsq=zsqinv*(hgtretc(k,ista)-hgtsng(jsta)) *              &
                            (hgtretc(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnretc(ista)-trnsng(jsta))*             &
                       (trnretc(ista)-trnsng(jsta))/(aglavg*aglavg)
                DO ivar=1,nvar
                  rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxret(ivar,k,ista)=oanxret(ivar,k,ista) +             &
                   (tabexp(indextab)*odifsng(ivar,jsta)*corsng(ivar,jsta))
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add contribution due to observation error variance
!
!-----------------------------------------------------------------------
!
        DO ivar=1,nvar
          oanxsng(ivar,jsta)=                                           &
                 oanxsng(ivar,jsta)+qobsng(ivar,jsta)*                  &
                 odifsng(ivar,jsta)*corsng(ivar,jsta)
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xsng(jsta))*(xradc(ista)-xsng(jsta))        &
                +(yradc(ista)-ysng(jsta))*(yradc(ista)-ysng(jsta))
            IF(dist < rlimsq) THEN
              DO k=1,nlevrad(ista)
                dzsq=zsqinv*(hgtradc(k,ista)-hgtsng(jsta)) *              &
                            (hgtradc(k,ista)-hgtsng(jsta))
                aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnradc(ista)-trnsng(jsta))*             &
                     (trnradc(ista)-trnsng(jsta))/(aglavg*aglavg)
!
                rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxrad(1,k,ista)=oanxrad(1,k,ista) +                     &
                  (tabexp(indextab)*odifsng(1,jsta)*corsng(1,jsta))
!
                rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxrad(2,k,ista)=oanxrad(2,k,ista) +                     &
                    (tabexp(indextab)*odifsng(2,jsta)*corsng(2,jsta))
!
                rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxrad(3,k,ista)=oanxrad(3,k,ista) +                     &
                    (tabexp(indextab)*odifsng(5,jsta)*corsng(5,jsta))
              END DO ! k
            END IF
          END IF
        END DO ! ista
      END IF ! source check
     END IF      ! isrcsng check
    END DO ! jsta
  END IF ! sngsw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
  IF(uasw > 0) THEN
    DO jsta=1,nobsua
      IF(iuseua(isrcua(jsta)) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xua(jsta))*(xsng(ista)-xua(jsta))            &
                +(ysng(ista)-yua(jsta))*(ysng(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtua(klev,jsta)) *             &
                            (hgtsng(ista)-hgtua(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnua(jsta))*               &
                           (trnsng(ista)-trnua(jsta))/(aglavg*aglavg)
            
                DO ivar=1,nvar
                  rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxsng(ivar,ista)=oanxsng(ivar,ista) +                 &
                    (tabexp(indextab)*odifua(ivar,klev,jsta)*             &
                         corua(ivar,klev,jsta))
                END DO ! ivar
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xua(jsta))*(xua(ista)-xua(jsta))              &
                +(yua(ista)-yua(jsta))*(yua(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnua(ista)-trnua(jsta))*                            &
                    (trnua(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtua(klev,jsta)) *          &
                              (hgtua(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)

                  DO ivar=1,nvar
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                      oanxua(ivar,k,ista)=oanxua(ivar,k,ista) +           &
                        (tabexp(indextab) * odifua(ivar,klev,jsta)*       &
                           corua(ivar,klev,jsta))
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add term due to observation error
!
!-----------------------------------------------------------------------
!
        DO klev=1,nlevsua(jsta)
          DO ivar=1,nvar
            oanxua(ivar,klev,jsta)=                                     &
                oanxua(ivar,klev,jsta)+qobsua(ivar,klev,jsta)*          &
                odifua(ivar,klev,jsta)*corua(ivar,klev,jsta)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to retrieval data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xua(jsta))*(xretc(ista)-xua(jsta))          &
                +(yretc(ista)-yua(jsta))*(yretc(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnua(jsta))*                          &
                    (trnretc(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtua(klev,jsta)) *        &
                              (hgtretc(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)

                  DO ivar=1,nvar
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                        oanxret(ivar,k,ista)=oanxret(ivar,k,ista) +       &
                         (tabexp(indextab)*odifua(ivar,klev,jsta)*        &
                           corua(ivar,klev,jsta))
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to radar data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xua(jsta))*(xradc(ista)-xua(jsta))          &
                +(yradc(ista)-yua(jsta))*(yradc(ista)-yua(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnua(jsta))*                          &
                    (trnradc(ista)-trnua(jsta))
              DO klev=1,nlevsua(jsta)
                agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtua(klev,jsta)) *        &
                              (hgtradc(k,ista)-hgtua(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(1,k,ista)=oanxrad(1,k,ista) +                  &
                    (tabexp(indextab)*odifua(1,klev,jsta)*               &
                         corua(1,klev,jsta))
!
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(2,k,ista)=oanxrad(2,k,ista) +                  &
                     (tabexp(indextab)*odifua(2,klev,jsta)*              &
                       corua(2,klev,jsta))
!
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(3,k,ista)=oanxrad(3,k,ista) +                  &
                       (tabexp(indextab)*odifua(5,klev,jsta)*            &
                           corua(5,klev,jsta))
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista

      END IF     ! source check
    END DO   ! ua jsta
  END IF   ! uasw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval data
!
!-----------------------------------------------------------------------
!
  IF(retsw > 0) THEN
    DO jsta=1,ncolret
      IF(iuseret(isrcret(iret(jsta))) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xretc(jsta))*(xsng(ista)-xretc(jsta))        &
                +(ysng(ista)-yretc(jsta))*(ysng(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtretc(klev,jsta)) *           &
                            (hgtsng(ista)-hgtretc(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnretc(jsta))*             &
                           (trnsng(ista)-trnretc(jsta))/(aglavg*aglavg)

                DO ivar=1,nvar
                  rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxsng(ivar,ista)=oanxsng(ivar,ista) +                 &
                     (tabexp(indextab)*odifret(ivar,klev,jsta)*           &
                         corret(ivar,klev,jsta))
                END DO ! ivar
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xretc(jsta))*(xua(ista)-xretc(jsta))          &
                +(yua(ista)-yretc(jsta))*(yua(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnua(jsta))*                          &
                    (trnretc(ista)-trnua(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtretc(klev,jsta))*         &
                              (hgtua(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)

                  DO ivar=1,nvar
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    oanxua(ivar,k,ista)=oanxua(ivar,k,ista) +             &
                       (tabexp(indextab)*odifret(ivar,klev,jsta)*         &
                           corret(ivar,klev,jsta))
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to retrieval data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xretc(jsta))*(xretc(ista)-xretc(jsta))      &
                +(yretc(ista)-yretc(jsta))*(yretc(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnretc(jsta))*                        &
                    (trnretc(ista)-trnretc(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtretc(klev,jsta)) *      &
                              (hgtretc(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)

                  DO ivar=1,nvar
                    rnorm=(dist*rngsqi(ivar)+dzsq+dtrnsq)*ftabinv
                    indextab=min((nint(rnorm)+1),ntabexp)
                    oanxret(ivar,k,ista)=oanxret(ivar,k,ista) +          &
                       (tabexp(indextab)*odifret(ivar,klev,jsta)*        &
                           corret(ivar,klev,jsta))
                  END DO ! ivar
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add term due to observation error
!
!-----------------------------------------------------------------------
!
        DO klev=1,nlevret(jsta)
          DO ivar=1,nvar
            oanxret(ivar,klev,jsta)=                                    &
                oanxret(ivar,klev,jsta)+qobsret(ivar,klev,jsta)*        &
                odifret(ivar,klev,jsta)*corret(ivar,klev,jsta)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xretc(jsta))*(xradc(ista)-xretc(jsta))      &
                +(yradc(ista)-yretc(jsta))*(yradc(ista)-yretc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnretc(jsta))*                        &
                    (trnradc(ista)-trnretc(jsta))
              DO klev=1,nlevret(jsta)
                agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtretc(klev,jsta)) *      &
                              (hgtradc(k,ista)-hgtretc(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(1,k,ista)=oanxrad(1,k,ista) +                  &
                    (tabexp(indextab)*odifret(1,klev,jsta)*              &
                         corret(1,klev,jsta))
!
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(2,k,ista)=oanxrad(2,k,ista) +                  &
                    (tabexp(indextab)*odifret(2,klev,jsta)*              &
                          corret(2,klev,jsta))
!
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(3,k,ista)=oanxrad(3,k,ista) +                  &
                    (tabexp(indextab)*odifret(5,klev,jsta)*              &
                         corret(5,klev,jsta))
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO !ista

      END IF ! source check
    END DO ! jsta
  END IF ! retsw check
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
  IF(radsw > 0) THEN
    DO jsta=1,ncolrad
      IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Add to single-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsng
         IF(isrcsng(ista) > 0) THEN
          IF(iusesng(isrcsng(ista)) > 0) THEN
            dist=(xsng(ista)-xradc(jsta))*(xsng(ista)-xradc(jsta))        &
                +(ysng(ista)-yradc(jsta))*(ysng(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                dzsq=zsqinv*(hgtsng(ista)-hgtradc(klev,jsta)) *           &
                            (hgtsng(ista)-hgtradc(klev,jsta))
                aglista=max(0.,(hgtsng(ista)-trnsng(ista)))
                aglavg=(0.5*(agljsta+aglista))+trnrcst
                dtrnsq=trnsqinv*(trnsng(ista)-trnradc(jsta))*             &
                       (trnsng(ista)-trnradc(jsta))/(aglavg*aglavg)
!
                rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxsng(1,ista)=oanxsng(1,ista) +                         &
                  (tabexp(indextab) * ABS(uazmrad(jsta))*                 &
                     odifrad(1,klev,jsta) * corrad(1,klev,jsta))
!
                rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxsng(2,ista)=oanxsng(2,ista) +                         &
                  (tabexp(indextab) * ABS(vazmrad(jsta))*                 &
                       odifrad(2,klev,jsta) * corrad(2,klev,jsta))
!
                rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                oanxsng(5,ista)=oanxsng(5,ista) +                         &
                   (tabexp(indextab) * odifrad(3,klev,jsta)*              &
                       corrad(3,klev,jsta))
              END DO ! klev
            END IF
          END IF
         END IF      ! isrcsng check
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to multiple-level data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,nobsua
          IF(iuseua(isrcua(ista)) > 0) THEN
            dist=(xua(ista)-xradc(jsta))*(xua(ista)-xradc(jsta))          &
                +(yua(ista)-yradc(jsta))*(yua(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnua(ista)-trnradc(jsta))*                          &
                    (trnua(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
                DO k=1,nlevsua(ista)
                  dzsq=zsqinv*(hgtua(k,ista)-hgtradc(klev,jsta))*         &
                       (hgtua(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtua(k,ista)-trnua(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxua(1,k,ista)=oanxua(1,k,ista) +                    &
                    (tabexp(indextab) * ABS(uazmrad(jsta))*              &
                       odifrad(1,klev,jsta) * corrad(1,klev,jsta))
!
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxua(2,k,ista)=oanxua(2,k,ista) +                    &
                    (tabexp(indextab) * ABS(vazmrad(jsta))*              &
                        odifrad(2,klev,jsta)*corrad(2,klev,jsta))
!
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxua(5,k,ista)=oanxua(5,k,ista) +                    &
                     (tabexp(indextab)*odifrad(3,klev,jsta)*             &
                         corrad(3,klev,jsta))
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
!-----------------------------------------------------------------------
!
!  Add to retrieval data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolret
          IF(iuseret(isrcret(iret(ista))) > 0) THEN
            dist=(xretc(ista)-xradc(jsta))*(xretc(ista)-xradc(jsta))      &
                +(yretc(ista)-yradc(jsta))*(yretc(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnretc(ista)-trnradc(jsta))*                        &
                    (trnretc(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))

                DO k=1,nlevret(ista)
                  dzsq=zsqinv*(hgtretc(k,ista)-hgtradc(klev,jsta)) *      &
                              (hgtretc(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtretc(k,ista)-trnretc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxret(1,k,ista)=oanxret(2,k,ista) +                  &
                    (tabexp(indextab) * ABS(uazmrad(jsta))*              &
                        odifrad(1,klev,jsta)*corrad(1,klev,jsta))
!
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxret(2,k,ista)=oanxret(2,k,ista) +                  &
                    (tabexp(indextab) * ABS(vazmrad(jsta))*              &
                        odifrad(2,klev,jsta)*corrad(2,klev,jsta))
!
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxret(5,k,ista)=oanxret(5,k,ista) +                  &
                    (tabexp(indextab)*                                   &
                        odifrad(3,klev,jsta)*corrad(3,klev,jsta))
                END DO
              END DO
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!
!  Add to radar data analyses
!
!-----------------------------------------------------------------------
!
        DO ista=1,ncolrad
          IF(iuserad(isrcrad(irad(ista))) > 0) THEN
            dist=(xradc(ista)-xradc(jsta))*(xradc(ista)-xradc(jsta))      &
                +(yradc(ista)-yradc(jsta))*(yradc(ista)-yradc(jsta))
            IF(dist < rlimsq) THEN
              dtrn2=(trnradc(ista)-trnradc(jsta))*                        &
                    (trnradc(ista)-trnradc(jsta))
              DO klev=1,nlevrad(jsta)
                agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))

                DO k=1,nlevrad(ista)
                  dzsq=zsqinv*(hgtradc(k,ista)-hgtradc(klev,jsta)) *      &
                              (hgtradc(k,ista)-hgtradc(klev,jsta))
                  aglista=max(0.,(hgtradc(k,ista)-trnradc(ista)))
                  aglavg=(0.5*(agljsta+aglista))+trnrcst
                  dtrnsq=trnsqinv*dtrn2/(aglavg*aglavg)
!
                  rnorm=(dist*rngsqi(1)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(1,k,ista)=oanxrad(1,k,ista) +                   &
                    (tabexp(indextab) * ABS(uazmrad(jsta))*               &
                        odifrad(1,klev,jsta) * corrad(1,klev,jsta))
!
                  rnorm=(dist*rngsqi(2)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(2,k,ista)=oanxrad(2,k,ista) +                   &
                    (tabexp(indextab) * ABS(vazmrad(jsta))*               &
                        odifrad(2,klev,jsta) * corrad(2,klev,jsta))
!
                  rnorm=(dist*rngsqi(5)+dzsq+dtrnsq)*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  oanxrad(3,k,ista)=oanxrad(3,k,ista) +                   &
                    (tabexp(indextab)*                                    &
                        odifrad(3,klev,jsta)*corrad(3,klev,jsta))
                END DO ! k
              END DO ! klev
            END IF
          END IF
        END DO ! ista
!
        DO klev=1,nlevrad(jsta)
          oanxrad(1,klev,jsta)=                                         &
                 oanxrad(1,klev,jsta)+qobsrad(1,klev,jsta)*             &
                 odifrad(1,klev,jsta)*corrad(1,klev,jsta)
          oanxrad(2,klev,jsta)=                                         &
                 oanxrad(2,klev,jsta)+qobsrad(2,klev,jsta)*             &
                 odifrad(2,klev,jsta)*corrad(2,klev,jsta)
          oanxrad(3,klev,jsta)=                                         &
                 oanxrad(3,klev,jsta)+qobsrad(3,klev,jsta)*             &
                 odifrad(3,klev,jsta)*corrad(3,klev,jsta)
        END DO
!
      END IF     ! source check
    END DO
  END IF   ! radsw check
  istatus=0
  RETURN
END SUBROUTINE brtobs3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BRTGRD3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE brtgrd3d(nx,ny,nz,nvar,nvarrad,                              &
           nzua,nzrdr,nzret,                                            &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,xs,ys,zs,icatg,xcor,     &
           odifsng,corsng,xsng,ysng,hgtsng,trnsng,                      &
           isrcsng,icatsng,nobsng,                                      &
           odifua,corua,xua,yua,hgtua,trnua,isrcua,nlevsua,nobsua,      &
           odifrad,corrad,xradc,yradc,hgtradc,trnradc,                  &
           uazmrad,vazmrad,irad,isrcrad,nlevrad,ncolrad,                &
           odifret,corret,xretc,yretc,hgtretc,trnretc,                  &
           iret,isrcret,nlevret,ncolret,                                &
           iusesng,iuseua,iuserad,iuseret,                              &
           sngsw,uasw,radsw,retsw,kpvrsq,rpass,                         &
           wlim,zrngsq,trnropt,trnrcst,trnrngsq,rngsqi,                 &
           anx,tem1,tem2,tem3,ibegin,iend,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Objectively analyze data on a regular grid specified
!  by x and y.
!  Bratseth scheme is employed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  July, 1995
!
!  MODIFICATION HISTORY:
!
!  Jan, 1996 (KB)
!  Added radar data and other improvements.
!
!  July, 1996 (KB)
!  Added sngsw and uasw to switch on/off use of single-level
!  and sounding data, respectively.
!
!  Jan 5, 1997 (KB)
!  New version to reorder summing to speed-up code.
!  Renamed to BRTGRD3D from BRTGRDZ3d and BRTGRDTH3d.
!  When called to analyze using delta-theta, potential
!  temperature data is passed into the height arrays.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Sizing arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nvarrad,nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
!
!-----------------------------------------------------------------------
!
!  Grid arguments
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  INTEGER :: icatg(nx,ny)
  REAL :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  REAL :: odifsng(nvar,mxsng)             ! variable to analyse
  REAL :: corsng(nvar,mxsng)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: trnsng(mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
!
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: corua(nvar,nzua,mxua)
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: trnua(mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: corrad(nvarrad,nzrdr,mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: trnradc(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: corret(nvar,nzret,mxcolret)
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: trnretc(mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Scratch space
!
!-----------------------------------------------------------------------
!
  REAL :: rngsqi(nvar)
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)
  INTEGER :: sngsw,uasw,radsw,retsw
  REAL :: kpvrsq(nvar)
  REAL :: rpass
  REAL :: wlim
  REAL :: zrngsq
  INTEGER :: trnropt
  REAL :: trnrcst
  REAL :: trnrngsq
!
!-----------------------------------------------------------------------
!
!  Output grid arguments
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  INTEGER :: ibegin(ny)
  INTEGER :: iend(ny)
  INTEGER :: istatus
!
  REAL :: ftabinv
  INTEGER :: ntabexp
  PARAMETER (ntabexp=5000)
  REAL :: tabexp(ntabexp)
  COMMON /ftabexp/ ftabinv
  COMMON /tablexp/ tabexp
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: rsq,rlim,rlimsq,dzsq,dtrnsq,zsqinv,trnsqinv,rnorm
  REAL :: agljsta,aglavg
  INTEGER :: i,j,k,jsta,klev,ivar,indextab,jbegin,jend

  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  rlimsq=0.
  DO ivar=1,nvar
    rsq=(kpvrsq(ivar)*rpass)
    rngsqi(ivar)=1./rsq
    rlimsq=AMAX1(rlimsq,rsq)
  END DO
  rlimsq=-rlimsq*ALOG(wlim)
  rlim=0.001*SQRT(rlimsq)
  IF (myproc == 0)                                                      &
    WRITE(6,'(/a,f10.2,a/)') ' Influence cutoff radius: ',rlim,' km.'
!
  zsqinv=1./zrngsq
!
  IF(trnropt > 0) THEN
    trnsqinv=1./trnrngsq
  ELSE
    trnsqinv=0.
  END IF


!-----------------------------------------------------------------------
!
!  Calculate Terrain and AGL heights
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1, nx
      tem1(i,j,2)=0.5*(zs(i,j,1)+zs(i,j,2))
    END DO
  END DO
  DO j=1,ny
    DO k=1, nz
      DO i=1, nx
        tem3(i,j,k)=max(0.,(zs(i,j,k)-tem1(i,j,2)))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  First sum contributions from single-level obs
!
!-----------------------------------------------------------------------
!
  IF(sngsw > 0) THEN
    DO jsta=1,nobsng
     IF(isrcsng(jsta) > 0) THEN
      IF(iusesng(isrcsng(jsta)) > 0) THEN
        agljsta=max(0.,(hgtsng(jsta)-trnsng(jsta)))
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,1)=(xs(i)-xsng(jsta))*(xs(i)-xsng(jsta))           &
                       +(ys(j)-ysng(jsta))*(ys(j)-ysng(jsta))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Find range in i,j of valid influence
!
!-----------------------------------------------------------------------
!
        jbegin=ny+1
        jend=1
        DO j=1,ny
          DO i=1,nx
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          ibegin(j)=i
          IF(i <= nx) THEN
            jbegin=MIN(jbegin,j)
            jend=MAX(jend,j)
          END IF
          DO i=nx,ibegin(j),-1
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          iend(j)=i
        END DO
!
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE (j,k,i)
        DO j=jbegin,jend
          DO k=1,nz-1
            DO i=ibegin(j),iend(j)
              dzsq=(zs(i,j,k)-hgtsng(jsta)) *                          &
                   (zs(i,j,k)-hgtsng(jsta)) * zsqinv
              aglavg=0.5*(agljsta+tem3(i,j,k))+trnrcst
              dtrnsq=trnsqinv*(trnsng(jsta)-tem1(i,j,2))*              &
                          (trnsng(jsta)-tem1(i,j,2))/(aglavg*aglavg)
              tem2(i,j,k)=dzsq+dtrnsq
            END DO
          END DO
        END DO
        DO ivar=1,nvar
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE (j,k,i,rnorm,indextab)
          DO j=jbegin,jend
            DO k=1,nz-1
              DO i=ibegin(j),iend(j)
                rnorm=(tem1(i,j,1)*rngsqi(ivar) + tem2(i,j,k))*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                anx(i,j,k,ivar)=anx(i,j,k,ivar) +                       &
                    (tabexp(indextab)*xcor(icatsng(jsta),icatg(i,j))    &
                    *odifsng(ivar,jsta)*corsng(ivar,jsta))
              END DO
            END DO
          END DO
        END DO
      END IF     ! source check
     END IF
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from multiple-level obs
!
!-----------------------------------------------------------------------
!
  IF(uasw > 0) THEN
    DO jsta=1,nobsua
      IF(iuseua(isrcua(jsta)) > 0) THEN
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,1)=(xs(i)-xua(jsta))*(xs(i)-xua(jsta))             &
                       +(ys(j)-yua(jsta))*(ys(j)-yua(jsta))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Find range in i,j of valid influence
!
!-----------------------------------------------------------------------
!
        jbegin=ny+1
        jend=1
        DO j=1,ny
          DO i=1,nx
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          ibegin(j)=i
          IF(i <= nx) THEN
            jbegin=MIN(jbegin,j)
            jend=MAX(jend,j)
          END IF
          DO i=nx,ibegin(j),-1
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          iend(j)=i
        END DO
!
        DO klev=1,nlevsua(jsta)
          agljsta=max(0.,(hgtua(klev,jsta)-trnua(jsta)))
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                dzsq=(zs(i,j,k)-hgtua(klev,jsta)) *                     &
                     (zs(i,j,k)-hgtua(klev,jsta)) * zsqinv
                aglavg=0.5*(agljsta+tem3(i,j,k))+trnrcst
                dtrnsq=trnsqinv*(trnua(jsta)-tem1(i,j,2))*              &
                          (trnua(jsta)-tem1(i,j,2))/(aglavg*aglavg)
                tem2(i,j,k)=dzsq+dtrnsq
              END DO
            END DO
          END DO
          DO ivar=1,nvar
            DO k=1,nz-1
              DO j=jbegin,jend
                DO i=ibegin(j),iend(j)
                  rnorm=(tem1(i,j,1)*rngsqi(ivar)+tem2(i,j,k))*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  anx(i,j,k,ivar)=anx(i,j,k,ivar) +                     &
                      (tabexp(indextab) *                               &
                       odifua(ivar,klev,jsta) *                         &
                       corua(ivar,klev,jsta))
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF     ! source check
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from radar obs
!
!-----------------------------------------------------------------------
!
  IF(radsw > 0) THEN
    DO jsta=1,ncolrad
      IF(iuserad(isrcrad(irad(jsta))) > 0) THEN
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,1)=(xs(i)-xradc(jsta))*(xs(i)-xradc(jsta))         &
                       +(ys(j)-yradc(jsta))*(ys(j)-yradc(jsta))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Find range in i,j of valid influence
!
!-----------------------------------------------------------------------
!
        jbegin=ny+1
        jend=1
        DO j=1,ny
          DO i=1,nx
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          ibegin(j)=i
          IF(i <= nx) THEN
            jbegin=MIN(jbegin,j)
            jend=MAX(jend,j)
          END IF
          DO i=nx,ibegin(j),-1
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          iend(j)=i
        END DO
!
        DO klev=1,nlevrad(jsta)
          agljsta=max(0.,(hgtradc(klev,jsta)-trnradc(jsta)))
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                dzsq=(zs(i,j,k)-hgtradc(klev,jsta)) *                   &
                     (zs(i,j,k)-hgtradc(klev,jsta)) * zsqinv
                aglavg=0.5*(agljsta+tem3(i,j,k))+trnrcst
                dtrnsq=trnsqinv*(trnradc(jsta)-tem1(i,j,2))*            &
                          (trnradc(jsta)-tem1(i,j,2))/(aglavg*aglavg)
                tem2(i,j,k)=dzsq+dtrnsq
              END DO
            END DO
          END DO
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                rnorm=(tem1(i,j,1)*rngsqi(1) + tem2(i,j,k))*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                anx(i,j,k,1)=anx(i,j,k,1) +                             &
                    (tabexp(indextab) *                                 &
                     ABS(uazmrad(jsta)) *                               &
                     odifrad(1,klev,jsta) *                             &
                     corrad(1,klev,jsta))
              END DO
            END DO
          END DO
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                rnorm=(tem1(i,j,1)*rngsqi(2) + tem2(i,j,k))*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                anx(i,j,k,2)=anx(i,j,k,2) +                             &
                       (tabexp(indextab) *                              &
                         ABS(vazmrad(jsta)) *                           &
                         odifrad(2,klev,jsta) *                         &
                         corrad(2,klev,jsta))
              END DO
            END DO
          END DO
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                rnorm=(tem1(i,j,1)*rngsqi(5) + tem2(i,j,k))*ftabinv
                indextab=min((nint(rnorm)+1),ntabexp)
                anx(i,j,k,5)=anx(i,j,k,5) +                             &
                       (tabexp(indextab) *                              &
                         odifrad(3,klev,jsta) *                         &
                         corrad(3,klev,jsta))
              END DO
            END DO
          END DO
        END DO
      END IF     ! source check
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Add contributions from retrieval data
!
!-----------------------------------------------------------------------
!
  IF(retsw > 0) THEN
    DO jsta=1,ncolret
      IF(iuseret(isrcret(iret(jsta))) > 0) THEN
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,1)=(xs(i)-xretc(jsta))*(xs(i)-xretc(jsta))         &
                       +(ys(j)-yretc(jsta))*(ys(j)-yretc(jsta))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!
!  Find range in i,j of valid influence
!
!-----------------------------------------------------------------------
!
        jbegin=ny+1
        jend=1
        DO j=1,ny
          DO i=1,nx
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          ibegin(j)=i
          IF(i <= nx) THEN
            jbegin=MIN(jbegin,j)
            jend=MAX(jend,j)
          END IF
          DO i=nx,ibegin(j),-1
            IF(tem1(i,j,1) < rlimsq) EXIT
          END DO
          iend(j)=i
        END DO
!
!-----------------------------------------------------------------------
!
!  Apply retrievel corrections
!
!-----------------------------------------------------------------------
!
        DO klev=1,nlevret(jsta)
          agljsta=max(0.,(hgtretc(klev,jsta)-trnretc(jsta)))
          DO k=1,nz-1
            DO j=jbegin,jend
              DO i=ibegin(j),iend(j)
                dzsq=(zs(i,j,k)-hgtretc(klev,jsta)) *                   &
                     (zs(i,j,k)-hgtretc(klev,jsta)) * zsqinv
                aglavg=0.5*(agljsta+tem3(i,j,k))+trnrcst
                dtrnsq=trnsqinv*(trnretc(jsta)-tem1(i,j,2))*            &
                          (trnretc(jsta)-tem1(i,j,2))/(aglavg*aglavg)
                tem2(i,j,k)=dzsq+dtrnsq
              END DO
            END DO
          END DO
          DO ivar=1,nvar
            DO j=jbegin,jend
              DO k=1,nz-1
                DO i=ibegin(j),iend(j)
                  rnorm=(tem1(i,j,1)*rngsqi(ivar)+tem2(i,j,k))*ftabinv
                  indextab=min((nint(rnorm)+1),ntabexp)
                  anx(i,j,k,ivar)=anx(i,j,k,ivar) +                     &
                           (tabexp(indextab) *                          &
                             odifret(ivar,klev,jsta) *                  &
                             corret(ivar,klev,jsta))
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF     ! source check
    END DO
  END IF
  istatus=0
  RETURN
END SUBROUTINE brtgrd3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE BARQC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE barqc(maxsta,nvar,nsrc,nprev,nsta,ipres,iptmp,iqv,           &
           obs,obstime,xsta,ysta,isrc,qual,                             &
           knt,wgtsum,zsum,sqsum,                                       &
           range,wlim,klim,qclim,varnam,snam,                           &
           iqclist,iqcsave,iqcflag)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Routine to objectively quality control surface data using a
!  Barnes analysis of nearby stations.
!
!  If an observations is found to be bad its quality variable
!  is set to the input flag value, iqcflag
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  September, 1989
!
!  MODIFICATION HISTORY:
!
!  Jul, 1996 (KB)
!  Version for ARPS analysis
!  Keith Brewster, CAPS
!
!  Aug, 1997 (KB)
!  Modifications and clarification of code for QC of analysis
!  variables instead of raw observed variables, before.  Related
!  to changes in data error specification in other routines.
!
!  May 2004 (KB)
!  Corrected checking of specific humidity, including calc of qvsat
!  and conversion of limits expressed in percent RH.  Added indices
!  of pres, potential temp and qv to argument list.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Observation Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxsta,nvar,nsrc
  INTEGER :: nprev,nsta
  INTEGER :: ipres,iptmp,iqv
!
  REAL :: obs(nvar,maxsta)
  INTEGER :: obstime(maxsta)
  REAL :: xsta(maxsta)
  REAL :: ysta(maxsta)
  INTEGER :: isrc(maxsta)
!
!-----------------------------------------------------------------------
!
!  Scratch space
!
!-----------------------------------------------------------------------
!
  INTEGER :: knt(nvar)
  REAL :: wgtsum(nvar)
  REAL :: zsum(nvar)
  REAL :: sqsum(nvar)
!
!-----------------------------------------------------------------------
!
!  File unit numbers
!
!-----------------------------------------------------------------------
!
  INTEGER :: iqclist,iqcsave
!
!-----------------------------------------------------------------------
!
!  Analysis specification arguments
!
!-----------------------------------------------------------------------
!
  REAL :: range     ! e-folding range km**2 of barnes weight
  REAL :: wlim      ! limit of weight to set max range
  INTEGER :: klim   ! minimum # of stations to influence grid pt
  REAL :: qclim(nvar,nsrc)   ! QC cutoffs
!
!-----------------------------------------------------------------------
!
!  Input/Output quality indicator
!
!-----------------------------------------------------------------------
!
  INTEGER :: qual(nvar,maxsta)
!
!-----------------------------------------------------------------------
!
!  Diagnostics and other argument variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: varnam(nvar)
  CHARACTER (LEN=5) :: snam(maxsta)
  INTEGER :: iqcflag
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: rlimsq,dist,wgt,dob,sqanx,chklim
  REAL :: pres,ptmp,temp,qvsat,qvlim

  INTEGER :: jsta,ksta,ivar
  LOGICAL :: listit,saveit
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
!  Set printing logicals.  Guaranteed MPI correct.  Only process 0
!  can have non-zero values.
!
!-----------------------------------------------------------------------
!
  listit=(iqclist > 0)
  saveit=(iqcsave > 0)
!
!-----------------------------------------------------------------------
!
!  Based on the minimum weight to consider, set
!  a maximum range
!
!-----------------------------------------------------------------------
!
  rlimsq=-range*ALOG(wlim)
!
!-----------------------------------------------------------------------
!
!   Uses Barnes weighting function.  see variable wgt.
!
!   Ksta is location that is being examined
!
!-----------------------------------------------------------------------
!
  PRINT *, ' BARQC: nsta= ',nsta
  DO ksta=(nprev+1),(nprev+nsta)
    IF(isrc(ksta) > 0) THEN
      DO ivar=1,nvar
        zsum(ivar)=0.
        sqsum(ivar)=0.
        wgtsum(ivar)=0.
        knt(ivar)=0
      END DO
      DO jsta=1,(nprev+nsta)
        IF(isrc(jsta) > 0) THEN
          IF(jsta /= ksta) THEN
            dist=(xsta(ksta)-xsta(jsta))*(xsta(ksta)-xsta(jsta))          &
                +(ysta(ksta)-ysta(jsta))*(ysta(ksta)-ysta(jsta))
            IF(dist < rlimsq) THEN
              wgt=EXP(-dist/range)
              DO ivar=1,nvar
                IF(qual(ivar,jsta) > 0) THEN
                  knt(ivar)=knt(ivar)+1
                  wgtsum(ivar)=wgtsum(ivar)+wgt
                  zsum(ivar)=zsum(ivar)+wgt*obs(ivar,jsta)
                  sqsum(ivar)=sqsum(ivar)+                                &
                       wgt*obs(ivar,jsta)*obs(ivar,jsta)
                END IF
              END DO
            END IF
          END IF  ! dont use current station
        END IF
      END DO
!
!  Normalization of weighted sums
!
      DO ivar=1,nvar
        IF(knt(ivar) > 0) THEN
          sqsum(ivar)=sqsum(ivar)-                                          &
                  (zsum(ivar)*zsum(ivar)/wgtsum(ivar))
          zsum(ivar)=zsum(ivar)/wgtsum(ivar)
        END IF
      END DO
!
!  Find appropriate temp and pres to compute qvsat for qv treshold
!
      IF( knt(ipres) >= klim ) THEN
        pres=zsum(ipres)
      ELSE IF (qual(ipres,ksta) > 0 ) THEN
        pres=obs(ipres,ksta)
      ELSE
        pres=p0
      END IF

      IF( knt(iptmp) >= klim ) THEN
        ptmp=zsum(iptmp)
      ELSE IF ( qual(iptmp,ksta) > 0 ) THEN
        ptmp=obs(iptmp,ksta)
      ELSE
        ptmp=293.0
      END IF
!
      temp=ptmp*((pres/p0)**rddcp)
      qvsat=f_qvsat(pres,temp)
!
!  Check all variables obs-analysis vs qc threshold.
!
      DO ivar=1,nvar
        IF(qual(ivar,ksta) > 0) THEN
          IF(knt(ivar) >= klim) THEN
            sqanx=MAX(sqsum(ivar),0.)
            sqanx=3.0*SQRT(sqanx/wgtsum(ivar))
            dob=obs(ivar,ksta)-zsum(ivar)

            IF(ivar /= iqv) THEN
              chklim=AMAX1(sqanx,qclim(ivar,isrc(ksta)))
            ELSE
              qvlim=qclim(iqv,isrc(ksta))*qvsat
              chklim=AMAX1(sqanx,qvlim)
            END IF

            IF(ABS(dob) > chklim) THEN
!
!  No checks for my proc here.  Any processor that has a complaint about an
!  ob that it handles needs to report it.
!
              PRINT *, ' QC flagging data...sta=',snam(ksta)
              PRINT *, '     var, ob =',                                &
                          varnam(ivar),obs(ivar,ksta)
              PRINT *, '     delta ob =',dob
              PRINT *, '  3.0*stdv,qclimit=',sqanx,                     &
                                   qclim(ivar,isrc(ksta))
!
!  "listit" and "saveit" are guaranteed to be always zero for all processors
!  except processor 0.
!
              IF(listit) WRITE(iqclist,810)                             &
                  snam(ksta),obstime(ksta),ivar,varnam(ivar),           &
                  obs(ivar,ksta),dob,                                   &
                  sqanx,qclim(ivar,isrc(ksta))
              810           FORMAT(2X,a5,i5,2X,i3,1X,a6,f11.2,f11.2,f11.2,f11.2)
              qual(ivar,ksta)=iqcflag
            END IF
            zsum(ivar)=dob
          ELSE
!          print *, ' Not enuf data to check ivar= ',
!    :                   ivar,' for ',snam(ksta)
            zsum(ivar)=-999.
          END IF
        ELSE
          zsum(ivar)=-999.
        END IF  ! data missing, dont bother checking
      END DO
      IF(saveit) WRITE(iqcsave,820) snam(ksta),obstime(ksta),           &
                (obs(ivar,ksta),ivar=1,nvar)
      820   FORMAT(2X,a5,i5,2X,15(f10.2))
    END IF     ! isrc check
  END DO
  RETURN
END SUBROUTINE barqc

!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION TEXP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION texp(rnorm)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A faster exponential function, via a lookup table.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Keith Brewster, CAPS
!  December, 1995
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
  REAL :: texp,setexp
  REAL :: rnorm
!
  REAL :: mxrnrm,ftabinv
  INTEGER :: ntabexp
  PARAMETER (ntabexp=5000)
  REAL :: tabexp(ntabexp)
  COMMON /ftabexp/ ftabinv
  COMMON /tablexp/ tabexp
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,indextab
  REAL :: ftab
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  indextab=nint(-(rnorm*ftabinv)) + 1
  indextab=MAX(indextab,1)
  indextab=MIN(indextab,ntabexp)
  texp=tabexp(indextab)
  RETURN
!
  ENTRY setexp(mxrnrm)
  ftab=mxrnrm/FLOAT(ntabexp)
  ftabinv=1./ftab
  DO i=1,ntabexp
    tabexp(i)=EXP((-(i-1)*ftab))
  END DO
!
! Force first element to one and last element to zero
!
  tabexp(1)=1.
  tabexp(ntabexp)=0.
  setexp=ftabinv
  RETURN
  END FUNCTION texp
