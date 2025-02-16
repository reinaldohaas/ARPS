!     ##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SOLAR_POS                   ######
!######             ARPS Data Analysis System                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solar_pos(rlat,rlon,i4time,alt,dec,hrangle)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine calculates the solar declination,
!  solar altitude angle, and the solar hour angle according to
!  the time and the latitude and longitude.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!          Based on the LAPS cloud analysis code by Steve Albers,
!          January,  1994.
!
!  MODIFICATION HISTORY:
!
!  04/22/96  J. Zhang
!            Modified for ADAS analysis.
!
!
!-----------------------------------------------------------------------
!
! Argument      I/O     Type                    Description
! --------      ---     ----    ----------------------------------------
! RLAT           I      R*4     Latitude (degrees)
! RLNG           I      R*4     Longitude (degrees)
! I4TIME         I       I      Time in absolute seconds from 00:00 UTC
!                            1/1/1960
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!  INPUT:
!
  REAL :: rlat,rlon            ! lat and lon of the data point.
  INTEGER :: i4time            ! the data time in seconds from 1/1/1977
!
!  OUTPUT:
!
  REAL :: alt      ! solar altitude angle for the data pt. at the time
  REAL :: dec      ! solar declination for the data pt. at the time
  REAL :: hrangle  ! hour angle for the data pt. at the time
!
!  LOCAL:
!
  INTEGER :: iyr,imon,idy,ihr,imin,isec,jd
  REAL :: pi,rpd
  PARAMETER (pi=3.1415926, rpd=pi/180.)
  REAL :: eqt,timeq,soldec,coszen
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
  CALL abss2ctim(i4time, iyr, imon, idy, ihr, imin, isec )
  CALL julday(iyr,imon,idy,jd)

  eqt = timeq(jd)/rpd            !Equation of Time (Degrees)
  dec = soldec(jd)/rpd           !Solar declination
  hrangle = (ihr-12)*15. + imin/4. + isec/240. + rlon + eqt

  coszen = SIN(rlat*rpd)*SIN(dec*rpd)                                   &
           + COS(rlat*rpd)*COS(dec*rpd)*COS(hrangle*rpd)
  alt = 90. - ACOS(coszen)/rpd

  IF(hrangle < -180.) hrangle = hrangle + 360.
  IF(hrangle > +180.) hrangle = hrangle - 360.

!    write(6,*)'jd,ihr,imin',jd,ihr,imin
!    write(6,*)'hrangle,dec,alt',hrangle,dec,alt

  RETURN
END SUBROUTINE solar_pos


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Thefollowing functions are used to calculate geometric solar
! parameters. These formulas are from Paltridge and Platt, 1976.
! They reference Spencer, 1971 for the solar declination and
! equation of time.
!
!
!  AUTHOR: J. Wakefield
!  01/28/1982, Original version
!
!  MODIFICATION HISTORY
!  04/22/1996 C Jian Zhang
!  Modified for ARPSDAS, added documents.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!

  FUNCTION radnorm(jd)
!
!-----------------------------------------------------------------------
!
! Purpose:
!
!  Calculate normalized earth-sun distance factor (R0/R)**2,
!  where JD is input Julian day number
!
!-----------------------------------------------------------------------

  dayang1=2.*3.14159265*(jd-1)/365.
  dayang2=2.*dayang1

  radnorm= 1.000110                                                     &
            +0.034221*COS(dayang1)+0.001280*SIN(dayang1)                &
            +0.000719*COS(dayang2)+0.000077*SIN(dayang2)

  RETURN
  END FUNCTION radnorm


!
!-----------------------------------------------------------------------
!

  FUNCTION soldec(jd)
!
!-----------------------------------------------------------------------
!
! Purpose:
!
!  Calculate solar declination angle (radians), JD is input
!  Julian day number
!
!-----------------------------------------------------------------------
!
  dayang1=2.*3.14159265*(jd-1)/365.
  dayang2=2.*dayang1
  dayang3=3.*dayang1

  soldec=  0.006918                                                     &
            -0.399912*COS(dayang1)+0.070257*SIN(dayang1)                &
            -0.006758*COS(dayang2)+0.000907*SIN(dayang2)                &
            -0.002697*COS(dayang3)+0.001480*SIN(dayang3)

  RETURN
  END FUNCTION soldec
!
!
!
!-----------------------------------------------------------------------
!

  FUNCTION timeq(jd)
!
!-----------------------------------------------------------------------
!
! Purpose:
!
!  Equation of time (radians), JD is input Julian day number
!
!-----------------------------------------------------------------------
!
  dayang1=2.*3.14159265*(jd-1)/365.
  dayang2=2.*dayang1

  timeq=   0.000075                                                     &
            +0.001868*COS(dayang1)-0.032077*SIN(dayang1)                &
            -0.014615*COS(dayang2)-0.040849*SIN(dayang2)

  RETURN
  END FUNCTION timeq
