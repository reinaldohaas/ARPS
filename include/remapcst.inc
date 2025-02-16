  INTEGER, PARAMETER :: maxlistf=10
  CHARACTER(LEN=4) :: radname
  CHARACTER(LEN=256) :: radfname
  CHARACTER(LEN=256) :: listfil(maxlistf)
  CHARACTER(LEN=80)  :: ref_time
  CHARACTER(LEN=80)  :: refvarname
  CHARACTER(LEN=80)  :: velvarname
  CHARACTER(LEN=80)  :: rhvvarname
  CHARACTER(LEN=80)  :: spwvarname
  CHARACTER(LEN=80)  :: snrvarname
  CHARACTER(LEN=80)  :: zdrvarname
  CHARACTER(LEN=80)  :: zdpvarname
  CHARACTER(LEN=80)  :: zvvvarname
  CHARACTER(LEN=80)  :: phivarname
  CHARACTER(LEN=80)  :: kdpvarname

  REAL :: vadradius,vadwidth
  REAL :: northazim,altoffset
  INTEGER :: nlistfil,radband,ncdcfile,rad98opt
  INTEGER :: qcopt,gclopt,velprocopt,dlpprocopt,clrvcpopt,tintrpopt
  INTEGER :: unfdiagopt,tradopt,fillref,medfilt,vadopt,vadminknt
  INTEGER :: nprocx_in,nprocy_in

  COMMON /radremap_opt/ radname,radfname,listfil,nlistfil,ref_time,    &
             refvarname,velvarname,rhvvarname,spwvarname,snrvarname,   &
             zdrvarname,zdpvarname,zvvvarname,phivarname,kdpvarname,   &
             radband,ncdcfile,rad98opt,qcopt,gclopt,                   &
             velprocopt,dlpprocopt,                                    &
             tintrpopt,clrvcpopt,unfdiagopt,tradopt,fillref,medfilt,   &
             northazim,altoffset,                                      &
             vadopt,vadradius,vadwidth,vadminknt,                      &
             nprocx_in,nprocy_in

  INTEGER :: iordref,iordvel,galcopt
  REAL :: refchek,refmiss
  REAL :: velchek,velmiss
  REAL :: rhvchek,rhvmiss
  REAL :: zdrchek,zdrmiss
  REAL :: kdpchek,kdpmiss
  REAL :: winszazim,winszrad
  REAL :: refmedl,velmedl,zdrmedl,kdpmedl,rhvmedl,refdazl,veldazl
  REAL :: spwthrrat,rhohvthr,snrthr,gcvrlim
  REAL :: dazim,rngmin,rngmax,envavgr,bmwidth

  COMMON /radremap_cst/ iordref,iordvel,                               &
             refchek,refmiss,velchek,velmiss,zdrchek,zdrmiss,          &
             kdpchek,kdpmiss,rhvchek,rhvmiss,winszazim,winszrad,       &
             refmedl,velmedl,zdrmedl,kdpmedl,rhvmedl,refdazl,veldazl,  &
             spwthrrat,rhohvthr,snrthr,gcvrlim,                        &
             dazim,rngmin,rngmax,envavgr,bmwidth

  INTEGER, PARAMETER :: maxfixed_nyqv = 10
  INTEGER :: nfixed_nyqv
  REAL :: fixed_nyqv(maxfixed_nyqv)
  REAL :: elevlim_fnyqv(maxfixed_nyqv)

  COMMON /radremap_nyqv/ nfixed_nyqv,fixed_nyqv,elevlim_fnyqv

  INTEGER :: fntimopt,dmpfmt,hdf4cmpr,dmpzero,dualpout
  REAL    :: outtime
  INTEGER :: ref2dopt,ref3dopt,vel2dopt,vel3dopt
  INTEGER :: zdr2dopt,zdr3dopt,kdp2dopt,kdp3dopt,rhv2dopt,rhv3dopt
  INTEGER :: wtsoloopt,wtsoloqcopt
  INTEGER :: wttiltopt,grdtiltver,wtgrdtiltopt
  COMMON /radremap_out/ fntimopt,dmpfmt,hdf4cmpr,                    &
          dmpzero,dualpout,outtime,                                  &
          ref2dopt,ref3dopt,vel2dopt,vel3dopt,                       &
          zdr2dopt,zdr3dopt,kdp2dopt,kdp3dopt,rhv2dopt,rhv3dopt,     &
          wttiltopt,wtgrdtiltopt,grdtiltver,wtsoloopt,wtsoloqcopt
