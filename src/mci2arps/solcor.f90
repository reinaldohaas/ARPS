SUBROUTINE solcor(britin,i4time,satlon,rlat,rlon,britout)
  IMPLICIT NONE
!
  INTEGER :: britin
  INTEGER :: i4time
  REAL :: satlon
  REAL :: rlat
  REAL :: rlon
  INTEGER :: britout
!
!  Satellite and sun position variables
!
  REAL :: latsat,lonsat,rangesat
  REAL :: pi,halfpi,dtr,rtd
  REAL :: solx,soly,solz
  REAL :: satx,saty,satz
  SAVE pi,halfpi,dtr,rtd
  SAVE solx,soly,solz
  SAVE satx,saty,satz
!
!  Constants used in photometric correction
!
  REAL :: s1,s2
  SAVE s1,s2
!
  REAL :: wfactor
  PARAMETER (wfactor=18.)
  REAL :: r01,rinf1,da1,r02,rinf2,da2,r31,r41,r32,r42
  PARAMETER (da1=12.,                                                   &
             r01=0.,                                                    &
             rinf1=71.,                                                 &
             da2=0.,                                                    &
             r02=80.,                                                   &
             rinf2=230.,                                                &
             r31=68.,                                                   &
             r41=220.,                                                  &
             r32=272.,                                                  &
             r42=880.)
  REAL :: r41mr31,r42mr32,a1,a2
  PARAMETER (r41mr31=r41-r31,                                           &
             r42mr32=r42-r32,                                           &
             a1=rinf1-r01,                                              &
             a2=rinf2-r02)
!
!  Physical constants
!
!  au      mean earth to sun distance (m)
!  eradius mean earth radius (m)
!  aunor   earth-to-sun distance normalized by earth radius
!
  REAL :: au,eradius,aunor
  PARAMETER (au = 1.4956E11,                                            &
             eradius = 6.371E06,                                        &
             aunor = au/eradius)
!
!  Misc local variables
!
  REAL :: alt,dec,hrangle
  REAL :: decrad,hrangrad,latrad,lonrad,altrad
  REAL :: pixx,pixy,pixz
  REAL :: rnor,delta,emiss,phase,ph2,phcor
  REAL :: r1,r2,rin11,rin12,rin21,rin22,rout1,rout2
!
!  Statement function
!
  REAL :: arctanh,x
  arctanh(x) = 0.5*LOG((1.+x)/(1.-x))
!
  CALL solar_pos(rlat,rlon,i4time,alt,dec,hrangle)
  IF(alt > 0.) THEN
!
!  Find position of pixel in cartesian coordinate
!  centered at center of earth, x axis to lat=0., lon=0.,
!  and zaxis to the north pole.
!  Earth's radius is used as unit length
!  so distance is unity for this calculation.
!
    latrad=dtr*rlat
    lonrad=dtr*rlon
    pixx=COS(lonrad)*COS(latrad)
    pixy=SIN(lonrad)*COS(latrad)
    pixz=SIN(latrad)
!
!  Find angle between vectors
!  satellite-to-pixel and center-of-earth-to-pixel
!
    CALL deltaang((satx-pixx),(saty-pixy),(satz-pixz),                  &
                   pixx,pixy,pixz,delta)

    emiss=halfpi-delta
!
!  Find angle between vectors
!  pixel-to-satellite and pixel-to-sun
!
    CALL deltaang((pixx-satx),(pixy-saty),(pixz-satz),                  &
                  (pixx-solx),(pixy-soly),(pixz-solz),delta)
    phase=AMAX1(AMIN1(halfpi,delta),0.)
!
!  Phase angle correction term
!
    altrad=dtr*alt
    ph2=COS(phase)
    ph2=ph2*ph2
    phcor=20.*(ph2*ph2*ph2)*SQRT(COS(emiss))*(SIN(altrad))**0.25
!
!  Stretching parameters
!
    r1=r01+a1*TANH(s1*dtr*(alt+da1))+phcor
    r2=r02+a2*TANH(s2*dtr*(alt+da2))+phcor
!
!  Apply stretching
!
    britout=nint(r1+(FLOAT(britin)-r1)*r41mr31/(r2-r1))
    britout=MIN(britout,255)
    britout=MAX(britout,0)
  ELSE
    britout=britin
  END IF
!
  RETURN
!
  ENTRY solr1r2(i4time,satlon,rlat,rlon,rout1,rout2)
  CALL solar_pos(rlat,rlon,i4time,alt,dec,hrangle)
  IF(alt > 0.) THEN
!
!  Find position of pixel in cartesian coordinate
!  centered at center of earth, x axis to lat=0., lon=0.,
!  and zaxis to the north pole.
!  Earth's radius is used as unit length
!  so distance is unity for this calculation.
!
    latrad=dtr*rlat
    lonrad=dtr*rlon
    pixx=COS(lonrad)*COS(latrad)
    pixy=SIN(lonrad)*COS(latrad)
    pixz=SIN(latrad)
!
!  Find angle between vectors
!  satellite-to-pixel and center-of-earth-to-pixel
!
    CALL deltaang((satx-pixx),(saty-pixy),(satz-pixz),                  &
                     pixx,pixy,pixz,delta)

    emiss=halfpi-delta
!
!  Find angle between vectors
!  pixel-to-satellite and pixel-to-sun
!
    CALL deltaang((pixx-satx),(pixy-saty),(pixz-satz),                  &
                  (pixx-solx),(pixy-soly),(pixz-solz),delta)
    phase=AMAX1(AMIN1(halfpi,delta),0.)
!
!  Phase angle correction term
!
    altrad=dtr*alt
    ph2=COS(phase)
    ph2=ph2*ph2
    phcor=20.*(ph2*ph2*ph2)*SQRT(COS(emiss))*(SIN(altrad))**0.25
!
!  Stretching parameters
!
    rout1=r01+a1*TANH(s1*dtr*(alt+da1))+phcor
    rout2=r02+a2*TANH(s2*dtr*(alt+da2))+phcor
  ELSE
    rout1=r31
    rout2=r41
  END IF
  RETURN

  ENTRY solrsc1(britin,rin11,rin21,britout)
  britout=nint(r31+(FLOAT(britin)-rin11)*r41mr31/(rin21-rin11))
  britout=MIN(britout,255)
  britout=MAX(britout,0)
  RETURN

  ENTRY solrsc2(britin,rin12,rin22,britout)
  britout=                                                              &
      nint(r32+(FLOAT(britin)-4*rin12)*r42mr32/(4*(rin22-rin12)))
  britout=MIN(britout,1023)
  britout=MAX(britout,0)
  RETURN
!
  ENTRY solcorset(i4time,latsat,lonsat,rangesat)
!
  pi=4.*ATAN(1.)
  dtr=pi/180.
  rtd=180./pi
  halfpi=0.5*pi
  s1=arctanh((r31-r01)/a1)/(dtr*(58.+da1))
  s2=arctanh((r41-r02)/a2)/(dtr*(58.+da2))
  PRINT *, ' Solcorset s1= ',s1,'  s2= ',s2
!
  CALL solar_pos(0.,0.,i4time,alt,dec,hrangle)
  decrad=dtr*dec
  hrangrad=dtr*hrangle
  solx=COS(-hrangrad)*COS(decrad)*aunor
  soly=SIN(-hrangrad)*COS(decrad)*aunor
  solz=SIN(decrad)*aunor
!
  latrad=dtr*latsat
  lonrad=dtr*lonsat
  rnor=rangesat/eradius
  satx=COS(lonrad)*COS(latrad)*rnor
  saty=SIN(lonrad)*COS(latrad)*rnor
  satz=SIN(latrad)*rnor

  RETURN
END SUBROUTINE solcor
!

SUBROUTINE deltaang(x1,y1,z1,x2,y2,z2,delta)
  IMPLICIT NONE
  REAL :: x1,y1,z1,x2,y2,z2,delta
!
  REAL :: vmag1,vmag2,vmag3,x3,y3,z3
  REAL :: dotp,sindel,cosdel
!
  vmag1=x1*x1+y1*y1+z1*z1
  vmag2=x2*x2+y2*y2+z2*z2
  x3=y1*z2 - z1*y2
  y3=z1*x2 - x1*z2
  z3=x1*y2 - y1*x2
  vmag3=x3*x3+y3*y3+z3*z3

  dotp=x1*x2+y1*y2+z1*z2

  sindel=vmag3/(vmag1*vmag2)
  cosdel=dotp/(vmag1*vmag2)

  delta=ATAN2(sindel,cosdel)

  RETURN
END SUBROUTINE deltaang
