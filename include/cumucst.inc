!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                     CUMUCST.INC                      ######
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
!  Include file 'cumucst.inc'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!  04/05/1995 (Yuhe Liu)
!  Documentation clean up.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nkp
  PARAMETER(nkp=80)

  INTEGER :: kcon,klcl,klfc,ketl,kct,igo,kmt

  REAL :: zmid,cdzmin,dzlow,dzhigh,plcl,tlcl,dzlcl

  COMMON/cpars /kcon,klcl,klfc,ketl,kct,igo,kmt,                        &
                zmid,cdzmin,dzlow,dzhigh,plcl,tlcl,dzlcl

  REAL :: cpr,alvl,aliv,aklv,akiv

  COMMON/const/   cpr,alvl,aliv,aklv,akiv

  REAL :: ucon(nkp),vcon(nkp),wcon(nkp),thtcon(nkp)                     &
                 ,qvcon(nkp),picon(nkp),tmpcon(nkp)                     &
                 ,dncon(nkp),zcon(nkp),zzcon(nkp)                       &
                 ,upe(nkp),vpe(nkp),wpe(nkp)                            &
                 ,ze(nkp),te(nkp),the(nkp),pe(nkp),rte(nkp)             &
                 ,pke(nkp),rhoe(nkp),thve(nkp),zc(nkp)                  &
                 ,rve(nkp),thee(nkp)                                    &
                 ,qvct1(nkp),qvct2(nkp),qvct3(nkp),qvct4(nkp)

  COMMON/cpenv / ucon,vcon,wcon,thtcon                                  &
                 ,qvcon,picon,tmpcon                                    &
                 ,dncon,zcon,zzcon                                      &
                 ,upe,vpe,wpe                                           &
                 ,ze,te,the,pe,rte                                      &
                 ,pke,rhoe,thve,zc                                      &
                 ,rve,thee                                              &
                 ,qvct1,qvct2,qvct3,qvct4

  REAL :: theu(nkp),rsu(nkp),thu(nkp),tu(nkp)                           &
                 ,thd(nkp),wtd(nkp),thcon(nkp),rtcon(nkp)

  COMMON/cpcon/   theu,rsu,thu,tu                                       &
                 ,thd,wtd,thcon,rtcon

  REAL :: vheat(nkp),vmois(nkp),vmdry(nkp)                              &
                 ,frcon(nkp),ftcon(nkp),cprecip,msflux                  &
                 ,tcon(nkp),rcon(nkp)

  COMMON/cptend/  vheat,vmois,vmdry                                     &
                 ,frcon,ftcon,cprecip,msflux                            &
                 ,tcon,rcon
