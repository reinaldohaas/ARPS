!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STATS                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE difstats(nx,ny,nz,                                           &
           nvar_anx,nvarradin,nvarrad,nzua,nzrdr,nzret,                 &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret,                          &
           xs,ys,zp,w,                                                  &
           xsng,ysng,hgtsng,thesng,                                     &
           obsng,odifsng,qobsng,qualsng,isrcsng,musesng,nobsng,         &
           xua,yua,hgtua,theua,                                         &
           obsua,odifua,qobsua,qualua,isrcua,museua,nlevsua,nobsua,     &
           elvrad,xradc,yradc,                                          &
           distrad,uazmrad,vazmrad,hgtradc,theradc,wradc,               &
           obsrad,odifrad,qobsrad,qualrad,                              &
           irad,isrcrad,muserad,nlevrad,ncolrad,                        &
           xretc,yretc,hgtretc,theretc,                                 &
           obsret,odifret,qobsret,qualret,                              &
           iret,isrcret,museret,nlevret,ncolret,                        &
           srcsng,srcua,srcrad,srcret,                                  &
           refmos,refgrid,                                              &
           knt,bias,rms,                                                &
           kntsngt,biassngt,rmssngt,kntuat,biasuat,rmsuat,              &
           kntrett,biasrett,rmsrett,kntradt,biasradt,rmsradt,           &
           kntsng,biassng,rmssng,kntua,biasua,rmsua,kntret,             &
           biasret,rmsret,kntrad,biasrad,rmsrad,                        &
           oanxsng,oanxua,oanxrad,oanxret,                              &
           tem1d,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE: Calculate summary tatistics from observation differences.
!
!  AUTHOR: Keith Brewster and Nazir Said, CAPS
!
!  INPUT: 
!  
!    nx         Number of grid points in the x-direction (east/west)
!    ny         Number of grid points in the y-direction (north/south)
!    nz         Number of grid points in the vertical
!    nvar_anx       Number of analysis variables
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
!    npass      Number of iterations
!    iwstat     Status indicator for writing statistics
!
!    xsng       x location of single-level data
!    ysng       y location of single-level data
!    hgtsng     elevation of single-level data  
!    thesng     theta (potential temperature) of single-level data
!    
!    obsng      single-level observations
!    odifsng    difference between single-level obs and analysis
!    qobsng     normalized observation error
!    qualsng    single-level data quality indicator
!    isrcsng    index of single-level data source
!    nobsng     number of single-level observations
!  
!    xua        x location of upper air data
!    yua        y location of upper air data
!    hgtua      elevation of upper air data
!    theua      theta (potential temperature) of upper air data
!
!    obsua      upper air observations
!    odifua     difference between upper air obs and analysis
!    qobsua     normalized observation error
!    qualua     upper air data quality indicator
!    isrcua     index of upper air data source
!    nlevsua    number of levels of data for each upper air location
!    nobsua     number of upper air observations
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
!    obsrad   radar observations
!    oanxrad  analysis (first guess) value at radar data location
!    odifrad  difference between radar observation and analysis
!    qobsrad  normalized observation error
!    qualrad  radar data quality indicator
!    ncolrad  number of radar columns read-in
!    istatus  status indicator
!    
!    latret   latitude of retrieval radar  (degrees N)
!    lonret   longitude of retrieval radar (degrees E)
!    xretc    x location of retrieval column
!    yretc    y location of retrieval column
!    iret     retrieval number
!    nlevret  number of levels of retrieval data in each column
!    hgtretc  height (m MSL) of retrieval observations
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
!
!-----------------------------------------------------------------------
!
!  Input Sizing Arguments
!    
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,lfnkey
  INTEGER :: nvar_anx,nvarradin,nvarrad,nvar
  INTEGER :: nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret,ncat
!
!-----------------------------------------------------------------------
!
!  Input gridded data
!
!-----------------------------------------------------------------------
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: w(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Input Observation Arguments 
!
!-----------------------------------------------------------------------
!
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: thesng(mxsng)
  REAL :: obsng(nvar_anx,mxsng)
  REAL :: odifsng(nvar_anx,mxsng)
  REAL :: qobsng(nvar_anx,mxsng)
  INTEGER :: qualsng(nvar_anx,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: musesng(0:nsrc_sng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
!
  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: theua(nzua,mxua)
  REAL :: obsua(nvar_anx,nzua,mxua)
  REAL :: odifua(nvar_anx,nzua,mxua)
  REAL :: qobsua(nvar_anx,nzua,mxua)
  INTEGER :: qualua(nvar_anx,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: museua(0:nsrc_ua)
  INTEGER :: nobsua
!
  REAL :: elvrad(mxrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: distrad(mxcolrad)
  REAL :: uazmrad(mxcolrad)   
  REAL :: vazmrad(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: theradc(nzrdr,mxcolrad)
  REAL :: wradc(nzrdr)
  REAL :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: qobsrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: muserad(0:nsrc_rad)
  INTEGER :: ncolrad
!
  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: theretc(nzret,mxcolret)
  REAL :: obsret(nvar_anx,nzret,mxcolret) 
  REAL :: odifret(nvar_anx,nzret,mxcolret)
  REAL :: qobsret(nvar_anx,nzret,mxcolret)
  INTEGER :: qualret(nvar_anx,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: museret(0:nsrc_ret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Input Analysis Control Variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcsng(nsrc_sng)
  CHARACTER (LEN=8) :: srcua (nsrc_ua )
  CHARACTER (LEN=8) :: srcrad(nsrc_rad)
  CHARACTER (LEN=8) :: srcret(nsrc_ret)
!
!-----------------------------------------------------------------------
!
! Input reflectivity data
!
!-----------------------------------------------------------------------
!
  REAL :: refmos(nx,ny,nz)
  REAL :: refgrid(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Scratch Space
!
!-----------------------------------------------------------------------
!
  INTEGER :: knt(nvar_anx)
  INTEGER :: kntsngt(nvar_anx)
  INTEGER :: kntuat(nvar_anx)
  INTEGER :: kntrett(nvar_anx)
  INTEGER :: kntradt(nvar_anx)
  INTEGER :: kntsng(nvar_anx,nsrc_sng)
  INTEGER :: kntua(nvar_anx,nsrc_ua)
  INTEGER :: kntret(nvar_anx,nsrc_ret)
  INTEGER :: kntrad(nvar_anx,nsrc_rad)
  REAL :: bias(nvar_anx)
  REAL :: rms(nvar_anx)
  REAL :: biassngt(nvar_anx)
  REAL :: rmssngt(nvar_anx)  
  REAL :: biasuat(nvar_anx)
  REAL :: rmsuat(nvar_anx)  
  REAL :: biasrett(nvar_anx)
  REAL :: rmsrett(nvar_anx)
  REAL :: biasradt(nvar_anx)
  REAL :: rmsradt(nvar_anx)
  REAL :: biassng(nvar_anx,nsrc_sng)
  REAL :: rmssng(nvar_anx,nsrc_sng)  
  REAL :: biasua(nvar_anx,nsrc_ua)
  REAL :: rmsua(nvar_anx,nsrc_ua)
  REAL :: biasret(nvar_anx,nsrc_ret)
  REAL :: rmsret(nvar_anx,nsrc_ret)
  REAL :: biasrad(nvar_anx,nsrc_rad)
  REAL :: rmsrad(nvar_anx,nsrc_rad)        
!
!-----------------------------------------------------------------------
! 
!  Output Variables at Observation Locations
!
!-----------------------------------------------------------------------
!
  REAL :: oanxsng(nvar_anx,mxsng)
  REAL :: oanxua(nvar_anx,nzua,mxua)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxret(nvar_anx,nzret,mxcolret)
!
  REAL :: tem1d(nz)
!
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(a,i4)') 'Computing orig difference stats'

  CALL verifst3d(nvar_anx,                                              &
           nzua,nzret,mxsng,mxua,mxret,mxcolret,                        &
           obsng,odifsng,oanxsng,isrcsng,musesng,qualsng,nobsng,        &
           obsua,odifua,oanxua,isrcua,museua,qualua,nlevsua,nobsua,     &
           obsret,odifret,oanxret,iret,isrcret,museret,                 &
           qualret,nlevret,ncolret,nsrc_sng,nsrc_ua,nsrc_ret,           &
           knt,bias,rms,kntsngt,biassngt,rmssngt,kntuat,biasuat,        &
           rmsuat,kntrett,biasrett,rmsrett,kntsng,biassng,rmssng,       &
           kntua,biasua,rmsua,kntret,biasret,rmsret)
!
  CALL verifstra(nx,ny,nz,nvar_anx,nvarrad,nvarradin,                   &
           nzrdr,mxrad,mxcolrad,                                        &
           xs,ys,zp,w,                                                  &
           elvrad,distrad,xradc,yradc,hgtradc,                          &
           uazmrad,vazmrad,wradc,qualrad,irad,isrcrad,                  &
           obsrad,odifrad,oanxrad,                                      &
           nlevrad,ncolrad,nsrc_rad,                                    &
           refmos,refgrid,                                              &
           kntradt,biasradt,rmsradt,kntrad,biasrad,rmsrad,tem1d)
  RETURN
END SUBROUTINE difstats

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VERIFST3d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE verifst3d(nvar_anx,                                          &
           nzua,nzret,mxsng,mxua,mxret,mxcolret,                        &
           obsng,odifsng,oanxsng,isrcsng,musesng,qualsng,nobsng,        &
           obsua,odifua,oanxua,isrcua,museua,qualua,nlevsua,nobsua,     &
           obsret,odifret,oanxret,iret,isrcret,museret,                 &
           qualret,nlevret,ncolret,nsrc_sng,nsrc_ua,nsrc_ret,           &
           knt,bias,rms,kntsngt,biassngt,rmssngt,kntuat,biasuat,        &
           rmsuat,kntrett,biasrett,rmsrett,kntsng,biassng,rmssng,       &
           kntua,biasua,rmsua,kntret,biasret,rmsret)

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
!  Dec, 2001 (KB & Nazir Said)
!  Updated
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nvar_anx,nz,nvar
  INTEGER :: nzua,nzret,mxsng,mxua,mxret,mxcolret
  INTEGER :: nsrc_sng,nsrc_ua,nsrc_ret
!
  REAL :: obsng(nvar_anx,mxsng)
  REAL :: odifsng(nvar_anx,mxsng)
  REAL :: oanxsng(nvar_anx,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: musesng(0:nsrc_sng)
  INTEGER :: qualsng(nvar_anx,mxsng)
  INTEGER :: nobsng
!
  REAL :: obsua(nvar_anx,nzua,mxua)
  REAL :: odifua(nvar_anx,nzua,mxua)
  REAL :: oanxua(nvar_anx,nzua,mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: museua(0:nsrc_ua)
  INTEGER :: qualua(nvar_anx,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: nobsua
!
  REAL :: obsret(nvar_anx,nzret,mxcolret)
  REAL :: odifret(nvar_anx,nzret,mxcolret)
  REAL :: oanxret(nvar_anx,nzret,mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(mxret)
  INTEGER :: museret(0:nsrc_ret) 
  INTEGER :: qualret(nvar_anx,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: ncolret
!
  INTEGER :: knt(nvar_anx)
  INTEGER :: kntsngt(nvar_anx)
  INTEGER :: kntuat(nvar_anx)
  INTEGER :: kntrett(nvar_anx)
  INTEGER :: kntsng(nvar_anx,nsrc_sng)
  INTEGER :: kntua(nvar_anx,nsrc_ua)
  INTEGER :: kntret(nvar_anx,nsrc_ret)
  REAL :: bias(nvar_anx)
  REAL :: rms(nvar_anx)
  REAL :: biassngt(nvar_anx)
  REAL :: biasuat(nvar_anx)
  REAL :: biasrett(nvar_anx)
  REAL :: rmssngt(nvar_anx)
  REAL :: rmsuat(nvar_anx)
  REAL :: rmsrett(nvar_anx)
  REAL :: biassng(nvar_anx,nsrc_sng)
  REAL :: rmssng(nvar_anx,nsrc_sng)
  REAL :: biasua(nvar_anx,nsrc_ua)
  REAL :: biasret(nvar_anx,nsrc_ret)
  REAL :: rmsua(nvar_anx,nsrc_ua)
  REAL :: rmsret(nvar_anx,nsrc_ret)

!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ivar,jsta,klev,jsrc
  REAL :: flknt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 DO ivar=1,nvar_anx
    knt(ivar) =0
    bias(ivar)=0.
    rms(ivar)=0.
    biassngt(ivar)=0.
    biasuat(ivar)=0.
    biasrett(ivar)=0.
    rmssngt(ivar)=0.
    rmsuat(ivar)=0.
    rmsrett(ivar)=0.
    kntsngt(ivar)=0
    kntuat(ivar)=0
    kntrett(ivar)=0
    DO jsrc=1,nsrc_sng
      kntsng(ivar,jsrc)=0
    END DO
    DO jsrc=1,nsrc_sng
      biassng(ivar,jsrc)=0.
      rmssng(ivar,jsrc)=0.
    END DO
    DO jsrc=1,nsrc_ua 
      kntua(ivar,jsrc)=0
      biasua(ivar,jsrc)=0.
      rmsua(ivar,jsrc)=0.
    END DO
    DO jsrc=1,nsrc_ret
      kntret(ivar,jsrc)=0
      biasret(ivar,jsrc)=0.
      rmsret(ivar,jsrc)=0.
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Find differences and form sums at single-level sites
!
!-----------------------------------------------------------------------
!
  DO jsta=1,nobsng
    jsrc=isrcsng(jsta)
    IF(jsrc /= 0) THEN
      IF(musesng(jsrc) > 0 ) THEN
        DO ivar=1,nvar_anx
          IF(qualsng(ivar,jsta) > 0) THEN
            odifsng(ivar,jsta)=(obsng(ivar,jsta)-oanxsng(ivar,jsta))
            knt(ivar)=knt(ivar)+1
            bias(ivar)=bias(ivar)-odifsng(ivar,jsta)
            rms(ivar)=rms(ivar)+                                        &
                       odifsng(ivar,jsta)*odifsng(ivar,jsta)
            kntsngt(ivar)=kntsngt(ivar)+1
            biassngt(ivar)=biassngt(ivar)-odifsng(ivar,jsta)
            rmssngt(ivar)=rmssngt(ivar)+                                &
                     (odifsng(ivar,jsta)*odifsng(ivar,jsta))
            kntsng(ivar,jsrc)=kntsng(ivar,jsrc)+1
            biassng(ivar,jsrc)=biassng(ivar,jsrc)-odifsng(ivar,jsta)
            rmssng(ivar,jsrc)= rmssng(ivar,jsrc)+                       &
                     (odifsng(ivar,jsta)*odifsng(ivar,jsta))
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
    jsrc=isrcua(jsta)
    IF(jsrc /= 0) THEN
      IF(museua(jsrc) > 0) THEN
        DO klev=1,nlevsua(jsta)
          DO ivar=1,nvar_anx
            IF(qualua(ivar,klev,jsta) > 0) THEN
              odifua(ivar,klev,jsta)=                                     &
                  obsua(ivar,klev,jsta)-oanxua(ivar,klev,jsta)
              knt(ivar)=knt(ivar)+1
              bias(ivar)=bias(ivar)-odifua(ivar,klev,jsta)
              rms(ivar)=rms(ivar)+                                        &
                odifua(ivar,klev,jsta)*odifua(ivar,klev,jsta)
              kntuat(ivar)=kntuat(ivar)+1
              biasuat(ivar)=biasuat(ivar)-odifua(ivar,klev,jsta)
              rmsuat(ivar)=rmsuat(ivar)+                                  &
                (odifua(ivar,klev,jsta)*odifua(ivar,klev,jsta))
              kntua(ivar,jsrc)=kntua(ivar,jsrc)+1
              biasua(ivar,jsrc)=biasua(ivar,jsrc)-odifua(ivar,klev,jsta)
              rmsua(ivar,jsrc)=rmsua(ivar,jsrc)+                          &
                (odifua(ivar,klev,jsta)*odifua(ivar,klev,jsta))
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
    jsrc=isrcret(iret(jsta))
    IF(jsrc > 0) THEN
      IF(museret(jsrc) > 0) THEN
        DO klev=1,nlevret(jsta)
          DO ivar=1,nvar_anx
            IF(qualret(ivar,klev,jsta) > 0) THEN
              odifret(ivar,klev,jsta)=                                    &
                  obsret(ivar,klev,jsta)-oanxret(ivar,klev,jsta)
              knt(ivar)=knt(ivar)+1
              bias(ivar)=bias(ivar)-odifret(ivar,klev,jsta)
              rms(ivar)=rms(ivar)+                                        &
                  odifret(ivar,klev,jsta)*odifret(ivar,klev,jsta)
              kntrett(ivar)=kntrett(ivar)+1
              biasrett(ivar)=biasrett(ivar)-odifret(ivar,klev,jsta)
              rmsrett(ivar)=rmsrett(ivar)+                                &
                 (odifret(ivar,klev,jsta)*odifret(ivar,klev,jsta))
              kntret(ivar,jsrc)=kntret(ivar,iret(jsta))+1
              biasret(ivar,jsrc)=biasret(ivar,jsrc)-odifret(ivar,klev,jsta)
              rmsret(ivar,jsrc)=rmsret(ivar,jsrc)+                        &
                 (odifret(ivar,klev,jsta)*odifret(ivar,klev,jsta))
            END IF
          END DO
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
  DO ivar=1,nvar_anx
    IF(knt(ivar) > 0) THEN
      flknt=FLOAT(knt(ivar))
      bias(ivar)=bias(ivar)/flknt
      rms(ivar)=SQRT(rms(ivar)/flknt)
    ELSE
      bias(ivar)=0.
      rms(ivar)=0.
    END IF
    IF(kntsngt(ivar) > 0) THEN
      flknt=FLOAT(kntsngt(ivar))
      biassngt(ivar)=biassngt(ivar)/flknt
      rmssngt(ivar)=SQRT(rmssngt(ivar)/flknt)
    ELSE
      biassngt(ivar)=0.
      rmssngt(ivar)=0.
    END IF
    IF(kntuat(ivar) > 0) THEN
      flknt=FLOAT(kntuat(ivar))
      biasuat(ivar)=biasuat(ivar)/flknt
      rmsuat(ivar)=SQRT(rmsuat(ivar)/flknt)
    ELSE
      biasuat(ivar)=0.
      rmsuat(ivar)=0.
    END IF
    IF(kntrett(ivar) > 0) THEN
      flknt=FLOAT(kntrett(ivar))
      biasrett(ivar)=biasrett(ivar)/flknt
      rmsrett(ivar)=SQRT(rmsrett(ivar)/flknt)
    ELSE
      biasrett(ivar)=0.
      rmsrett(ivar)=0.
    END IF
    
    DO jsrc =1,nsrc_sng
      IF(kntsng(ivar,jsrc) > 0 ) THEN
        flknt=FLOAT(kntsng(ivar,jsrc))
        biassng(ivar,jsrc)=biassng(ivar,jsrc)/flknt
        rmssng(ivar,jsrc)=SQRT(rmssng(ivar,jsrc)/flknt)
      ELSE
        biassng(ivar,jsrc)=0.
        rmssng(ivar,jsrc)=0.
      END IF
    END DO
    DO jsrc =1,nsrc_ua
      IF(kntua(ivar,jsrc) > 0 ) THEN
        flknt=FLOAT(kntua(ivar,jsrc))
        biasua(ivar,jsrc)=biasua(ivar,jsrc)/flknt
        rmsua(ivar,jsrc)=SQRT(rmsua(ivar,jsrc)/flknt)
      ELSE
        biasua(ivar,jsrc)=0.
        rmsua(ivar,jsrc)=0.
      END IF
    END DO    
    DO jsrc =1,nsrc_ret
      IF(kntret(ivar,jsrc) > 0 ) THEN
        flknt=FLOAT(kntret(ivar,jsrc))
        biasret(ivar,jsrc)=biasret(ivar,jsrc)/flknt
        rmsret(ivar,jsrc)=SQRT(rmsret(ivar,jsrc)/flknt)
      ELSE
        biasret(ivar,jsrc)=0.
        rmsret(ivar,jsrc)=0.
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE verifst3d

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VERIFSTRA                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE verifstra(nx,ny,nz,nvar_anx,nvarrad,nvarradin,               &
           nzrdr,mxrad,mxcolrad,                                        &
           xs,ys,zp,w,                                                  &
           elvrad,distrad,xradc,yradc,hgtradc,                          &
           uazmrad,vazmrad,wradc,qualrad,irad,isrcrad,                  &
           obsrad,odifrad,oanxrad,                                      &
           nlevrad,ncolrad,nsrc_rad,                                    &
           refmos,refgrid,                                              &
           kntradt,biasradt,rmsradt,kntrad,biasrad,rmsrad,hgt1d)
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
!  August, 2003 (KB)
!  Updated for statistical calculations
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nvar_anx,nvarrad,nvarradin,nzrdr,mxrad,mxcolrad,nsrc_rad
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zp(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: elvrad(mxrad)
  REAL :: distrad(mxcolrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  REAL :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: wradc(nzrdr)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: ncolrad
!
  REAL :: refmos(nx,ny,nz)
  REAL :: refgrid(nx,ny,nz)
!
  INTEGER :: kntradt(nvar_anx)
  REAL :: biasradt(nvar_anx)
  REAL :: rmsradt(nvar_anx)
  INTEGER :: kntrad(nvar_anx,nsrc_rad)
  REAL :: biasrad(nvar_anx,nsrc_rad)
  REAL :: rmsrad(nvar_anx,nsrc_rad)

  REAL :: hgt1d(nz)
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: icol,ivar,kr,jrad,jsrc
  REAL :: pi,dtr,dz,eleva,range,flknt
  REAL :: vr,vt,vrdiff,avrdif,dhdr,dsdr,dsdrinv,refdiff
  REAL :: vr_miss,refobs,refgrd
!
!-----------------------------------------------------------------------
!
!  Radar compariason parameters
!
!  refmiss   Lower limit of observed reflectivity
!  refmin    Lower limit of meaningful radar reflectivity
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refmiss=-21.
  REAL, PARAMETER :: refmin=10.
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
  DO ivar=1,nvar_anx
    kntradt(ivar)=0.
    biasradt(ivar)=0.
    rmsradt(ivar)=0.
    DO jsrc=1,nsrc_rad
     kntrad(ivar,jsrc)=0
     biasrad(ivar,jsrc)=0.
     rmsrad(ivar,jsrc)=0.
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Find radial velocity differences at radar data points
!  Add to sums
!
!-----------------------------------------------------------------------
!
  print *, ' Inside radar stats, ncolrad=',ncolrad
  DO icol=1,ncolrad
    IF(irad(icol) > 0) THEN
      jrad=irad(icol)
      jsrc=isrcrad(jrad)
!
      IF(jsrc > 0) THEN

        CALL wintrp(nx,ny,nz,nzrdr,nlevrad, xs,ys,zp,w,                    &
             xradc(icol),yradc(icol),hgtradc(1,icol),wradc,hgt1d)
!
        DO kr=1,nlevrad(icol)
          IF(qualrad(1,kr,icol) > 0 .AND. qualrad(2,kr,icol) > 0  .AND.    &
             hgtradc(kr,icol) > hgt1d(2) .AND.                             &
             hgtradc(kr,icol) < hgt1d(nz-1) ) THEN
!
            dz=hgtradc(kr,icol)-elvrad(jrad)
            CALL beamelv(dz,distrad(icol),eleva,range)
            CALL dhdrange(eleva,range,dhdr)
            dsdr=SQRT(MAX(0.,(1.-dhdr*dhdr)))
            IF(dsdr /= 0.) THEN
              dsdrinv=1./dsdr
            ELSE
              dsdrinv=0.
            END IF
            CALL vterm(obsrad(1,kr,icol),hgtradc(kr,icol),vt)
!
! Obtain w from grid, which was not previously saved.
!
            vr=((uazmrad(icol)*oanxrad(1,kr,icol) +                        &
                 vazmrad(icol)*oanxrad(2,kr,icol)) * dsdr) +               &
                 ((wradc(kr)-vt) * dhdr )
            vrdiff=vr-obsrad(2,kr,icol)
            avrdif=ABS(vrdiff)
            IF(icol < 50) print *, ' vr, vrobs, vrdiff=',vr,obsrad(2,kr,icol),vrdiff
            IF(icol < 50) print *, ' u,v,w=',oanxrad(1,kr,icol),oanxrad(2,kr,icol),wradc(kr)
            odifrad(1,kr,icol)=vrdiff
            kntradt(1)=kntradt(1)+1
            biasradt(1)=biasradt(1)+vrdiff
            rmsradt(1)=rmsradt(1)+(vrdiff*vrdiff)
            kntrad(1,jsrc)=kntrad(1,jsrc)+1
            biasrad(1,jsrc)=biasrad(1,jsrc)+vrdiff
            rmsrad(1,jsrc)=rmsrad(1,jsrc)+(vrdiff*vrdiff)
          END IF
        END DO
      END IF
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate errors in reflectivity.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-1
    DO j=2,ny-2
      DO i=2,nx-2
        IF(refmos(i,j,k) > refmiss) THEN
          refobs=max(refmos(i,j,k),refmin)
          refgrd=max(refgrid(i,j,k),refmin)
          refdiff=refgrd-refobs
          IF(abs(refdiff) > 60.) THEN
            print *, ' refgrd, refobs: ',refgrd,refobs
            print *, '         orig  : ',refmos(i,j,k),refgrid(i,j,k)
          END IF
          kntradt(2)=kntradt(2)+1
          biasradt(2)=biasradt(2)+refdiff
          rmsradt(2)=rmsradt(2)+(refdiff*refdiff)
        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate errors in composite reflectivity
!
!-----------------------------------------------------------------------
!
  DO j=2,ny-2
    DO i=2,nx-2
      IF(refmos(i,j,1) > refmiss) THEN
        refobs=max(refmos(i,j,1),refmin)
        refgrd=max(refgrid(i,j,1),refmin)
        refdiff=refgrd-refobs
        IF(abs(refdiff) > 60.) THEN
          print *, ' refgrd, refobs: ',refgrd,refobs
          print *, '         orig  : ',refmos(i,j,1),refgrid(i,j,1)
        END IF
        kntradt(3)=kntradt(3)+1
        biasradt(3)=biasradt(3)+refdiff
        rmsradt(3)=rmsradt(3)+(refdiff*refdiff)
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate stats from the sums formed above
!
!-----------------------------------------------------------------------
!
  DO ivar=1,nvar_anx
    IF(kntradt(ivar) > 0) THEN
      flknt=FLOAT(kntradt(ivar))
      biasradt(ivar)=biasradt(ivar)/flknt
      flknt=FLOAT(kntradt(ivar)-1)
      rmsradt(ivar)=SQRT(rmsradt(ivar)/flknt)
    ELSE
      biasradt(ivar)=0.
      rmsradt(ivar)=0.
    END IF
    DO jsrc=1,nsrc_rad 
      IF(kntrad(ivar,jsrc) > 0) THEN
       flknt=FLOAT(kntrad(ivar,jsrc))
       biasrad(ivar,jsrc)=biasrad(ivar,jsrc)/flknt
       rmsrad(ivar,jsrc)=SQRT(rmsrad(ivar,jsrc)/flknt)
      ELSE
       biasrad(ivar,jsrc)=0.
       rmsrad(ivar,jsrc)=0.
      END IF
    END DO
  END DO
!
  RETURN
END SUBROUTINE verifstra
!
SUBROUTINE wintrp(nx,ny,nz,nzrdr,nlevrad, xs,ys,zp,w,    &
             xradc,yradc,hgtradc,wradc,hgt1d)
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nzrdr
  INTEGER :: nlevrad
  REAL :: xs(nx)
  REAL :: ys(nx)
  REAL :: zp(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: xradc
  REAL :: yradc
  REAL :: hgtradc(nzrdr)
  REAL :: wradc(nzrdr)
  REAL :: hgt1d(nz)
!
! Misc local variables
!
  INTEGER :: i,j,k,kr,imid,jmid,kmid,iloc,jloc,kloc
  REAL :: wx,wy,wlo,whi
  REAL :: w11,w12,w21,w22
!
  imid=nx/2
  jmid=ny/2
  kmid=nz/2
  wradc=0.0
!
!-----------------------------------------------------------------------
!
!  Find column in grid
!
!-----------------------------------------------------------------------
!
  IF(xradc < xs(imid)) THEN
    DO i=imid,2,-1
      IF(xs(i) < xradc) EXIT
    END DO
    iloc=i
  ELSE
    DO i=imid,nx-1
      IF(xs(i) > xradc) EXIT
    END DO
    iloc=i-1
  END IF
  wx =  (xradc-xs(iloc))/(xs(iloc+1)-xs(iloc))
!
  IF(yradc < ys(jmid)) THEN
    DO j=jmid,2,-1
      IF(ys(j) < yradc) EXIT
    END DO
    jloc=j
  ELSE
    DO j=jmid,ny-1
      IF(ys(j) > yradc) EXIT
    END DO
    jloc=j-1
  END IF
  wy = (yradc-ys(jloc))/(ys(jloc+1)-ys(jloc))
!
!-----------------------------------------------------------------------
!
!  Determine bilinear interpolation weights
!
!-----------------------------------------------------------------------
!
  w22 = wx*wy
  w12 = (1.0 - wx) * wy
  w21 = wx * (1.0 - wy)
  w11 = (1.0 - wx) * (1.0 - wy)
!
!-----------------------------------------------------------------------
!
!  Interpolate grid heights in horizontal
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    hgt1d(k)=w11*zp(  iloc,  jloc,k) +                                  &
             w21*zp(iloc+1,  jloc,k) +                                  &
             w12*zp(  iloc,jloc+1,k) +                                  &
             w22*zp(iloc+1,jloc+1,k)
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop in height
!
!-----------------------------------------------------------------------
!
  DO kr=1,nlevrad
    IF(hgtradc(kr) >= hgt1d(2) .AND. hgtradc(kr) <= hgt1d(nz-1)) THEN
!
!-----------------------------------------------------------------------
!
!  Find z location in grid
!
!-----------------------------------------------------------------------
!
      IF(hgtradc(kr) < hgt1d(kmid)) THEN
        DO k=kmid,2,-1
          IF(hgt1d(k) < hgtradc(kr)) EXIT
        END DO
        kloc=k
      ELSE
        DO k=kmid,nz-1
          IF(hgt1d(k) > hgtradc(kr)) EXIT
        END DO
        kloc=k-1
      END IF
!
!-----------------------------------------------------------------------
!
!  Set z weights
!
!-----------------------------------------------------------------------
!
      whi = (hgtradc(kr)-hgt1d(kloc))/                                  &
            (hgt1d(kloc+1)-hgt1d(kloc))
      wlo = 1.0 - whi
      wradc(kr)=                                                        &
          wlo* (w11*w(  iloc,  jloc,  kloc) +                           &
                w21*w(iloc+1,  jloc,  kloc) +                           &
                w12*w(  iloc,jloc+1,  kloc) +                           &
                w22*w(iloc+1,jloc+1,  kloc))                            &
        + whi* (w11*w(  iloc,  jloc,kloc+1) +                           &
                w21*w(iloc+1,  jloc,kloc+1) +                           &
                w12*w(  iloc,jloc+1,kloc+1) +                           &
                w22*w(iloc+1,jloc+1,kloc+1))
    END IF
  END DO
  RETURN
END SUBROUTINE wintrp

SUBROUTINE vterm(obsref,hgt,vt)
  IMPLICIT NONE
  REAL :: obsref
  REAL :: hgt
  REAL :: vt

  REAL, PARAMETER :: zfrez=3000.
  REAL, PARAMETER :: zice=8000.
  REAL, PARAMETER :: rho0=1.2250
  REAL, PARAMETER :: h0=7000.
  REAL, PARAMETER :: denom=(1./(zice-zfrez))
  REAL, PARAMETER :: dbzmin=10.
  REAL, PARAMETER :: dbzmax=90.
!
  REAL :: refz,rhofact,s1,s2
!
  IF(obsref > dbzmin .AND. obsref < dbzmax ) THEN
    refz=10.**(0.1*obsref)
    rhofact=EXP(0.4*hgt/h0)
    IF(hgt < zfrez) THEN
      vt=2.6*(refz**0.107)*rhofact
    ELSE IF(hgt < zice) THEN
      s1=(zice-hgt)*denom
      s2=2.*(hgt-zfrez)*denom
      vt=s1*2.6*(refz**0.107)*rhofact + s2
    ELSE
      vt=2.0
    END IF
  ELSE 
    vt=0.
  END IF
  RETURN
END SUBROUTINE vterm
