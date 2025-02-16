!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_TIME_STRING            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_time_string  (curtime, time_string, tzabbr, tzone)

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  REAL,    INTENT(IN) :: curtime
  CHARACTER (LEN=25)  :: time_string
  INTEGER, INTENT(IN) :: tzone          ! time zone offset in hours
  CHARACTER (LEN=3), INTENT(IN) :: tzabbr

!-----------------------------------------------------------------------

  INTEGER :: abstsec
  INTEGER :: rjday, wday
  INTEGER :: ryear,rmonth, rday, rhour, rmin, rsec
  CHARACTER (LEN=3) :: smonth(12)
  CHARACTER (LEN=3) :: weekday(0:6)
  CHARACTER (LEN=2) :: chour, cmin

  DATA smonth/ 'Jan','Feb','Mar','Apr','May','Jun',                     &
               'Jul','Aug','Sep','Oct','Nov','Dec'/
  DATA weekday/'Sun','Mon','Tue','Wed','Thu','Fri','Sat'/

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL  ctim2abss( year,month,day,hour,minute,second, abstsec )

  !print*,year,month,day,hour,minute,second, abstsec
  !print*,'current time', curtime

  abstsec = abstsec + INT(curtime) + tzone*3600

  CALL abss2ctim( abstsec, ryear, rmonth, rday, rhour, rmin, rsec)
  !print*,abstsec, ryear, rmonth, rday, rhour, rmin, rsec

  CALL julday( ryear, rmonth, rday, rjday )
  CALL getwekday ( ryear,rmonth, rday, wday )

  WRITE(chour,'(i2.2)') rhour
  WRITE(cmin, '(i2.2)') rmin

  WRITE(time_string,'(8a,I2.2,3a,I4.4)') chour,':', cmin,' ',TRIM(tzabbr),' ', &
                    weekday(wday),' ',rday,' ',smonth(rmonth),' ',ryear

  RETURN
END SUBROUTINE get_time_string

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETWEKDAY                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getwekday ( year, month, day, wekday )

  INTEGER :: year, month, day, wekday
  INTEGER :: jday, tmp1
  INTEGER :: baseyear
  baseyear = 1960   ! have to be leap year

  CALL julday( year, month, day, jday )  ! current jday

  tmp1 = year-baseyear   ! total years

  lp = (year + 3 - baseyear) / 4         ! Number of leap days since 1960

  jday = jday + tmp1*365 + lp  ! total day from baseyear/1/1/

  IF( year == baseyear .AND. jday >= 60) jday = jday+1

  wekday = MOD(jday+4, 7)   !! 1960 1/1 is Friday

  RETURN
END SUBROUTINE getwekday

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_forecast_time          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_forecast_hms(timsnd, timhms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given a forecast time in seconds, convert it into a string for forecast
!  time in the hour/minute/second format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  06/28/2013
!  Rewritten from and replaced cvttim in src/arps/outlib3d.f90.
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    timsnd   Time in seconds
!
!  OUTPUT:
!
!    timhms    string contain time in hour:minute:second format
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  REAL, INTENT(IN)  :: timsnd  ! Time in seconds
  CHARACTER(LEN=9)  :: timhms  ! string contain time in
                               ! hour:minute:second format
  INTEGER :: h,m,s
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  h = INT(timsnd/3600.0 )
  m = INT((timsnd-h*3600.0)/60.0)
  s = nint(timsnd-h*3600.0-m*60.0)

  IF( s == 60) THEN
    m = m+1
    s = 0
  END IF

  WRITE(timhms,'(I0,a,i2.2,a,i2.2)') h,':',m,':',s

  RETURN
END SUBROUTINE get_forecast_hms