!
!##################################################################
!##################################################################
!######                                                      ######
!######          SUBROUTINE MEAN_ABSOULTE_ERROR              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE MEAN_ABSOLUTE_ERROR (nx, ny, nz, data, model, ilevel, &
                                jlevel, klevel)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculates the mean absolute error between the hypothetical data set
!  and model forecast and writes this to the "MAE.txt" file
!
!  The file "MAE.txt" will always contain the mean absolute error for 
!  the entire data set at the top of the file.  It will also contain
!  the mean absolute error based on the i, j, or k levels depending
!  upon which of those options are chosen in the "verif.input" file.
!
!  Formula:
!
!  MAE = (1/n) * SUM[k=1 to n] ABS(y(k) - o(k))
!
!  n is the number of forecast / observation pairs
!  y(k) is the kth of n forecasts
!  o(k) is the kth of n observations               
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Geoffrey Stano
!  03/01/2002 Initialized.
!
!  MODIFICATIONS:
!
!  04/04/2002 Jason J. Levit
!  Cosmetic modifications.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

       IMPLICIT none

       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios

       REAL :: maez(nz)
       REAL :: maey(ny)
       REAL :: maex(nx)
       REAL :: mae_all   !  variable for MAE of entire data set 
       REAL :: num
       REAL :: data(nx,ny,nz), model(nx,ny,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!  ##### Calculates mean absolute error for entire data set mae_all #####

       mae_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                num = 0.0
                num = model(i,j,k) - data(i,j,k)
                IF (num < 0) THEN
                   num = -1 * num  !  Performs Absolute Value
                END IF
                mae_all = mae_all + num
             END DO
          END DO
       END DO

       OPEN (UNIT = 11, FILE = "MAE.txt", IOSTAT = ios, &
       FORM = "FORMATTED")

       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF

       WRITE (11, FMT = 200)

       mae_all = mae_all / (nx * ny * nz)
       WRITE (11,'(f7.2)') mae_all

! ########## End of mae_all calculation ##########


! ##### Calculates i-level mean absolute error, meax ##### 

       IF (ilevel == 1) THEN
       DO i=1,nx
          maex(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                num = 0.0
                num = model(i,j,k) - data(i,j,k)
                IF (num < 0) THEN
                   num = -1 * num  !  Performs Absolute Value
                END IF
                maex(i) = maex(i) + num
             END DO
          END DO
       END DO
       
       WRITE (11, FMT = 201)

       DO i=1,nx   
          maex(i) = maex(i) / (ny * nz)
          WRITE (11,'(f7.2)') maex(i)
       END DO
       END IF

! ########## End of i-level mean absolute error ##########

! ##### Calculates j-level mean absolute error, maey #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          maey(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                num = 0.0
                num = model(i,j,k) - data(i,j,k)
                IF (num < 0) THEN
                   num = -1 * num  !  Performs Absolute Value
                END IF
                maey(j) = maey(j) + num
             END DO
          END DO
       END DO

       WRITE (11, FMT = 202)
       
       DO j=1,ny 
          maey(j) = maey(j) / (nx * nz)
          WRITE (11,'(f7.2)') maey(j)
       END DO
       END IF

! ######### End of j-level mean absolute error ##########

! ##### Calculates the k-level mean absolute error #####

       IF (klevel == 1) THEN
       DO k=1,nz
          maez(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                num = 0.0
                num = model(i,j,k) - data(i,j,k)
                IF (num < 0) THEN
                   num = -1 * num  !  Performs Absolute Value
                END IF
                maez(k) = maez(k) + num
             END DO
          END DO
       END DO

       WRITE (11, FMT = 203)
       
       DO k=1,nz 
          maez(k) = maez(k) / (nx * ny)
          WRITE (11,'(f7.2)') maez(k)
       END DO
       END IF

! ######### End of k-level mean absolute error #########

       CLOSE (UNIT = 11)


200    FORMAT ("Mean Absolute Error for all data"/)
201    FORMAT (/"Mean Absolute Error with i-level dependence"/)
202    FORMAT (/"Mean Absolute Error with j-level dependence"/)
203    FORMAT (/"Mean Absolute Error with k-level dependence"/)

       END SUBROUTINE Mean_Absolute_Error



!#######################################################################

       SUBROUTINE BRIER_SCORE (noe, fcstprob, obsevent)


!  Calculates the brier score between the hypothetical probability
!  forecast and the observed event and writes this to the "BrierS.txt" 
!  file
!
!  Formula:
!
!  BS = (1/n) * SUM[k=1 to n] (fp(k) - v(k))**2
!
!  n is the number of forecast / observation pairs 
!  fp is the forecasted probability for an event 
!  v is the verification term 
!  v = 1 if the event did occur
!  v = 0 if the event did not occur

!#######################################################################

       IMPLICIT none
          
       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
     
       REAL :: brier_s(noe)   !  brier score
       
       REAL :: fcstprob(noe), obsevent(noe)

       DO g=1,noe
          brier_s = brier_s + (fcstprob(g) - obsevent(g))**2
       END DO

       OPEN (UNIT = 22, FILE = "BrierS.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF

       WRITE (22, FMT = 300)

       brier_s(noe) = (brier_s(noe) / noe)
       WRITE (22, FMT = '(f4.2)') brier_s(noe)
   
       CLOSE(22)

300    FORMAT ("Brier Score"/)     


       END SUBROUTINE Brier_Score




!#######################################################################

       SUBROUTINE ROOT_MEAN_SQUARE (nx, ny, nz, data, model, ilevel, &
        jlevel, klevel,rms_all)


!  Calculates the root mean square error between the hypothetical data set
!  and the model forecast and writes this to the "RMS.txt" file
!
!  Formula:
! 
!  RMS = SQRT{ (1/n) * SUM[k=1 to n] (y(k) - o(k))**2}
!
!  n is the number of forecast / observation pairs 
!  y(k) is the kth of n forecasts
!  o(k) is the kth of n observations

!#######################################################################

       IMPLICIT none

       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num
       INTEGER :: rmscount

       REAL :: rms_all
       REAL :: rmsz(nz)
       REAL :: rmsy(ny)
       REAL :: rmsx(nx)
       REAL :: data(nx,ny,nz), model(nx,ny,nz)

! ##### Calulates the root mean square error for all data #####

       rms_all = 0.0
       rmscount=0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                IF (data(i,j,k).ne.-99.9) THEN
                  rms_all = rms_all + (model(i,j,k)-data(i,j,k))**2
                  rmscount=rmscount+1
                END IF
             END DO
          END DO
       END DO

       IF ( rmscount == 0) THEN
          rms_all = -99.9
          RETURN
       END IF

!       OPEN (UNIT = 12, FILE = "RMS.txt", IOSTAT = ios, &
!       FORM = "FORMATTED")
   
!      IF (ios /= 0) THEN
!         PRINT*, "Error opening file"
!      END IF
 
!       WRITE (12, FMT = 300)
!       print *, 'rmscount=',rmscount
!       rms_all = rms_all / (nx * ny * nz)
       rms_all = rms_all / real(rmscount)
       rms_all = SQRT(rms_all)
!       WRITE (12,'(f7.2)') rms_all

! ########## End all root mean square error ##########

! ##### Calculates the i-level root mean square error #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          rmsx(i) = 0.0
          DO j=1,ny
             DO k=1,nz  
                rmsx(i) = rmsx(i) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
!       WRITE (12, FMT = 301)

       DO i = 1,nx
          rmsx(i) = rmsx(i) / (ny * nz)
          rmsx(i) = SQRT(rmsx(i))
!          WRITE (12,'(f7.2)') rmsx(i)
       END DO
       END IF

! ########## End i-level root mean square error ##########

! ##### Calculates the j-level root mean square error #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          rmsy(j) = 0.0
          DO i=1,nx
             DO k=1,nz  
                rmsy(j) = rmsy(j) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
!       WRITE (12, FMT = 302)

       DO j = 1,ny
          rmsy(j) = rmsy(j) / (nx * nz)
          rmsy(j) = SQRT(rmsy(j))
!          WRITE (12,'(f7.2)') rmsy(j)
       END DO
       END IF

! ########## End j-level root mean square error ##########

! ##### Calculates the k-level root mean square error #####

       IF (klevel == 1) THEN
       DO k=1,nz
          rmsz(k) = 0.0
          DO j=1,ny
             DO i=1,nx  
                rmsz(k) = rmsz(k) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
!       WRITE (12, FMT = 303)

       DO k = 1,nz
          rmsz(k) = rmsz(k) / (nx * ny)
          rmsz(k) = SQRT(rmsz(k))
!          WRITE (12,'(f7.2)') rmsz(k)
       END DO
       END IF

! ########## End k-level root mean square error ##########
       
!       CLOSE (UNIT = 12)

300    FORMAT ("Root Mean Square Error for all data"/)
301    FORMAT (/"Root Mean Square Error with i-level dependence"/)
302    FORMAT (/"Root Mean Square Error with j-level dependence"/)
303    FORMAT (/"Root Mean Square Error with k-level dependence"/)


       END SUBROUTINE Root_Mean_Square



!#######################################################################

       SUBROUTINE MEAN_SQUARE_ERROR (nx, ny, nz, data, model, ilevel, &
         jlevel, klevel)


!  Calculates the mean square error between the hypothetical data set 
!  and the model forecast and writes this to the "MSE.txt" file
! 
!  Formula:
!
!  MSE = (1/n) * SUM[k=1 to n] (y(k) -o(k))**2
!
!  n is the number of forecast / observation pairs
!  y(k) is the kth of n forecasts
!  o(k) is the kth of n observations

!#######################################################################

       IMPLICIT none

       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num

       REAL :: mse_all
       REAL :: msez(nz)
       REAL :: msey(ny)
       REAL :: msex(nx)
       REAL :: data(nx,ny,nz), model(nx,ny,nz)
          
! ##### Calculates the mean square error for all data #####

       mse_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                mse_all = mse_all + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO

       OPEN (UNIT = 13, FILE = "MSE.txt", IOSTAT = ios,  &
       FORM = "FORMATTED")
   
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF

       WRITE (13, FMT = 400)
       mse_all = mse_all / (nx * ny * nz)
       WRITE (13,'(f7.2)') mse_all

! ########## End mean square error for all data ##########

! ##### Calculates i-level mean square error #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          msex(i) = 0.0
          DO j=1,ny   
             DO k=1,nz
                msex(i) = msex(i) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
       WRITE (13, FMT = 401)

       DO i = 1,nx
          msex(i) = msex(i) / (ny * nz)
          WRITE (13,'(f7.2)') msex(i)
       END DO
       END IF

! ######### End i-level mean square error ##########

! ##### Calculates j-level mean square error #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          msey(j) = 0.0
          DO k=1,nz   
             DO i=1,nx
                msey(j) = msey(j) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
       WRITE (13, FMT = 402)

       DO j = 1,ny
          msey(j) = msey(j) / (nx * nz)
          WRITE (13,'(f7.2)') msey(j)
       END DO
       END IF

! ########## End j-level mean square error ##########

! ##### Calculates k-level mean square error #####

       IF (klevel == 1) THEN
       DO k=1,nz
          msez(k) = 0.0
          DO j=1,ny   
             DO i=1,nx
                msez(k) = msez(k) + (model(i,j,k)-data(i,j,k))**2
             END DO
          END DO
       END DO
       
       WRITE (13, FMT = 403)

       DO k = 1,nz
          msez(k) = msez(k) / (nx * ny)
          WRITE (13,'(f7.2)') msez(k)
       END DO
       END IF

! ######### End k-level mean square error ##########
       
       CLOSE (UNIT = 13)

400    FORMAT ("Mean Square Error for all data"/)
401    FORMAT (/"Mean Square Error with i-level dependence"/)
402    FORMAT (/"Mean Square Error with j-level dependence"/)
403    FORMAT (/"Mean Square Error with k-level dependence"/)
       
          
       END SUBROUTINE Mean_Square_Error



!#######################################################################

       SUBROUTINE FORECAST_MEAN (nx, ny, nz, model, ilevel, jlevel, &
         klevel)


!  Calculates the mean of all of the forecast data and writes this 
!  to the "FcstMean.txt" file
!
!  Formula:
!
!  Forecast Mean = (Sum of all forecasts) / (Total number of forecasts)

!#######################################################################

       IMPLICIT none

       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num

       REAL :: fcstmean_all   !  mean of all forecasts
       REAL :: fcstmeanz(nz)
       REAL :: fcstmeany(ny)
       REAL :: fcstmeanx(nx)
       REAL :: model(nx,ny,nz)   !  hypothetical model data
       
! ##### Calculates the forecast mean for all data #####

      fcstmean_all = 0.0 

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                fcstmean_all = fcstmean_all + model(i,j,k)
             END DO
          END DO
       END DO

       fcstmean_all = fcstmean_all / (nx * ny * nz)

       OPEN (UNIT = 14, FILE = "FcstMean.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
       
       WRITE (14, FMT = 500)
       WRITE (14,'(f7.2)') fcstmean_all

! ############ End forecast mean for all data ##########

! ##### Calculates i-level forecast mean #####

       IF (ilevel == 1) THEN
       DO i=1,nx      
          fcstmeanx(i) = 0.0
          DO j=1,ny   
             DO k=1,nz
                fcstmeanx(i) = fcstmeanx(i) + model(i,j,k)
             END DO
          END DO
       END DO

       WRITE (14, FMT = 501)

       DO i=1,nx       
          fcstmeanx(i) = fcstmeanx(i) / (ny * nz)
          WRITE (14,'(f7.2)') fcstmeanx(i)
       END DO
       END IF

! ########## End i-level forecast mean ##########

! ##### Calculates j-level forecast mean #####

       IF (jlevel == 1) THEN
       DO j=1,ny      
          fcstmeany(j) = 0.0
          DO k=1,nz   
             DO i=1,nx
                fcstmeany(j) = fcstmeany(j) + model(i,j,k)
             END DO
          END DO
       END DO

       WRITE (14, FMT = 502)

       DO j=1,ny       
          fcstmeany(j) = fcstmeany(j) / (nx * nz)
          WRITE (14,'(f7.2)') fcstmeany(j)
       END DO
       END IF

! ########## End j-level forecast mean ##########

! ##### Calculates k-level forecast mean #####
       
       IF (klevel == 1) THEN
       DO k=1,nz      
          fcstmeanz(k) = 0.0
          DO j=1,ny   
             DO i=1,nx
                fcstmeanz(k) = fcstmeanz(k) + model(i,j,k)
             END DO
          END DO
       END DO

       WRITE (14, FMT = 503)

       DO k=1,nz       
          fcstmeanz(k) = fcstmeanz(k) / (nx * ny)
          WRITE (14,'(f7.2)') fcstmeanz(k)
       END DO
       END IF

! ########## End k-level forecast mean ##########

       CLOSE (UNIT = 14)

500    FORMAT ("Forecast mean for all data"/)
501    FORMAT (/"Forecast mean with i-level dependence"/)
502    FORMAT (/"Forecast mean with j-level dependence"/)
503    FORMAT (/"Forecast mean with k-level dependence"/)


       END SUBROUTINE Forecast_Mean




!#######################################################################

       SUBROUTINE OBSERVATION_MEAN (nx, ny, nz, data, ilevel, jlevel, &
        klevel)


!  Calculates the mean of all of the observation data and writes this
!  to the "ObsMean.txt" file
!
!  Formula:
! 
!  Observation Mean = (Sum of all obs) / (Total number of obs)

!#######################################################################
       
       IMPLICIT none
       
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num
          
       REAL :: obsmean_all   !  mean of all observations
       REAL :: obsmeanz(nz)
       REAL :: obsmeany(ny)
       REAL :: obsmeanx(nx)
       REAL :: data(nx,ny,nz)   !  hypothetical observation data
          
! ##### Calculates the mean for all the observations #####

       obsmean_all = 0.0
   
       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                obsmean_all = obsmean_all + data(i,j,k)
             END DO
          END DO   
       END DO      
       
       obsmean_all = obsmean_all / (nx * ny * nz)
       
       OPEN (UNIT = 15, FILE = "ObsMean.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
     
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
       
       WRITE (15, FMT = 600)
       WRITE (15,'(f7.2)') obsmean_all

! ######### End observation mean for all data ##########

! ##### Calculates i-level observation mean #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          obsmeanx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                obsmeanx(i) = obsmeanx(i) + data(i,j,k)
             END DO
          END DO
       END DO

       WRITE (15, FMT = 601)

       DO i=1,nx
          obsmeanx(i) = obsmeanx(i) / (ny * nz)
          WRITE (15,'(f7.2)') obsmeanx(i)
       END DO
       END IF

! ########## End i-level observation mean #########

! ##### Calculates j-level observation mean #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          obsmeany(j) = 0.0
          DO k=1,nz
             DO i=1,nx
                obsmeany(j) = obsmeany(j) + data(i,j,k)
             END DO
          END DO
       END DO

       WRITE (15, FMT = 602)

       DO j=1,ny
          obsmeany(j) = obsmeany(j) / (nx * nz)
          WRITE (15,'(f7.2)') obsmeany(j)
       END DO
       END IF

! ########## End j-level observation mean ##########

! ##### Calculates k-level observation mean #####

       IF (klevel == 1) THEN
       DO k=1,nz
          obsmeanz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                obsmeanz(k) = obsmeanz(k) + data(i,j,k)
             END DO
          END DO
       END DO

       WRITE (15, FMT = 603)

       DO k=1,nz
          obsmeanz(k) = obsmeanz(k) / (nx * ny)
          WRITE (15,'(f7.2)') obsmeanz(k)
       END DO
       END IF

! ########## End k-level observation mean ##########
       
       CLOSE (UNIT = 15)

600    FORMAT ("Observation mean for all data"/)
601    FORMAT (/"Observation mean with i-level dependence"/)
602    FORMAT (/"Observation mean with j-level dependence"/)
603    FORMAT (/"Observation mean with k-level dependence"/)
       
       
       END SUBROUTINE Observation_Mean
       
       


!#######################################################################

       SUBROUTINE FORECAST_VARIANCE (nx, ny, nz, model, ilevel, &
         jlevel, klevel)


!  Calculates the varience for the hypothetical forecast set and writes 
!  this to the "FcstVar.txt" file.
!
!  Formula:
!
!  Variance = SUM[ (x - x(bar))**2 / N]
!
!  x is the value of each individual forecast
!  x(bar) is the mean of the set of forecasts
!  N is the number of forecasts

!#######################################################################

       IMPLICIT none
       
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num

       REAL :: avg_all   !  holds average value
       REAL :: avgz(nz)
       REAL :: avgy(ny)
       REAL :: avgx(nx)
       REAL :: fcstvar_all
       REAL :: fcstvarz(nz)
       REAL :: fcstvary(ny)
       REAL :: fcstvarx(nx)
       REAL :: model(nx,ny,nz)

! ##### Calculates forecast variance for all data #####

       avg_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                avg_all = avg_all + model(i,j,k)
             END DO
          END DO
       END DO

       avg_all = avg_all / (nx * ny * nz)

       fcstvar_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                fcstvar_all = fcstvar_all + (model(i,j,k) - avg_all)**2
             END DO
          END DO
       END DO
       
       fcstvar_all = fcstvar_all / (nx*ny*nz)

       OPEN (UNIT = 17, FILE = "FcstVar.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
       
       WRITE (17, FMT = 700)
       WRITE (17,'(e13.2)') fcstvar_all

! ########## End forecast variance for all data ##########
          
! ##### Calculates i-level forecast variance #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          avgx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                avgx(i) = avgx(i) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO i=1,nx
          avgx(i) = avgx(i) / (ny * nz)
       END DO   
    
       DO i=1,nx
          fcstvarx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                fcstvarx(i) = fcstvarx(i) + (model(i,j,k) - avgx(i))**2
             END DO
          END DO
       END DO

       WRITE (17, FMT = 701)

       DO i=1,nx
          fcstvarx(i) = fcstvarx(i) / (ny*nz)
          WRITE (17,'(e13.2)') fcstvarx(i)
       END DO
       END IF

! ########## End i-level forecast variance ##########

! ##### Calculates j-level forecast variance #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          avgy(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                avgy(j) = avgy(j) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO j=1,ny
          avgy(j) = avgy(j) / (nx * nz)
       END DO
          
       DO j=1,ny
          fcstvary(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                fcstvary(j) = fcstvary(j) + (model(i,j,k) - avgy(j))**2
             END DO
          END DO
       END DO
       
       WRITE (17, FMT = 702)
       
       DO j=1,ny
          fcstvary(j) = fcstvary(j) / (nx * nz)
          WRITE (17,'(e13.2)') fcstvary(j)
       END DO
       END IF

! ########## End j-level forecast variance ##########

! ##### Calculates k-level forecast variance #####

       IF (klevel == 1) THEN
       DO k=1,nz
          avgz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                avgz(k) = avgz(k) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO k=1,nz
          avgz(k) = avgz(k) / (nx * ny)
       END DO
          
       DO k=1,nz
          fcstvarz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                fcstvarz(k) = fcstvarz(k) + (model(i,j,k) - avgz(k))**2
             END DO
          END DO
       END DO
       
       WRITE (17, FMT = 703)
       
       DO k=1,nz
          fcstvarz(k) = fcstvarz(k) / (nx * ny)
          WRITE (17,'(e13.2)') fcstvarz(k)
       END DO
       END IF

! ########## End k-level forecast variance ##########

       CLOSE (17)

700    FORMAT ("Forecast Variance with all data"/)
701    FORMAT (/"Forecast Variance with i-level dependence"/)
702    FORMAT (/"Forecast Variance with j-level dependence"/)
703    FORMAT (/"Forecast Variance with k-level dependence"/)


       END SUBROUTINE FORECAST_VARIANCE



!#######################################################################
       
       SUBROUTINE OBSERVATION_VARIANCE (nx, ny, nz, obs, ilevel, &
         jlevel, klevel)
       
          
!  Calculates the varience for the hypothetical forecast set and writes
!  this to the "FcstVar.txt" file.
!
!  Formula:
!
!  Variance = SUM[ (x - x(bar))**2 / N]
!
!  x is the value of each individual observation
!  x(bar) is the mean of the set of observations
!  N is the number of observations
             
!#######################################################################

       IMPLICIT none
       
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation   
       INTEGER :: i, j ,k, ios, num
          
       REAL :: avg_all   !  holds average value
       REAL :: avgz(nz)
       REAL :: avgy(ny)
       REAL :: avgx(nx)
       REAL :: obsvar_all
       REAL :: obsvarz(nz)
       REAL :: obsvary(ny)
       REAL :: obsvarx(nx)
       REAL :: obs(nx,ny,nz), model(nx,ny,nz)

! ##### Calculates observation variance for all data #####

       avg_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                avg_all = avg_all + obs(i,j,k)
             END DO
          END DO
       END DO

       avg_all = avg_all / (nx * ny * nz)

       obsvar_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                obsvar_all = obsvar_all + (obs(i,j,k) - avg_all)**2
             END DO
          END DO
       END DO
       
       obsvar_all = obsvar_all / (nx*ny*nz)

       OPEN (UNIT = 18, FILE = "ObsVar.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
       
       WRITE (18, FMT = 800)
       WRITE (18,'(e13.2)') obsvar_all

! ########## End observation variance for all data ##########
          
! ##### Calculates i-level observation variance #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          avgx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                avgx(i) = avgx(i) + obs(i,j,k)
             END DO
          END DO
       END DO
       
       DO i=1,nx
          avgx(i) = avgx(i) / (ny * nz)
       END DO   
    
       DO i=1,nx
          obsvarx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                obsvarx(i) = obsvarx(i) + (obs(i,j,k) - avgx(i))**2
             END DO
          END DO
       END DO

       WRITE (18, FMT = 801)

       DO i=1,nx
          obsvarx(i) = obsvarx(i) / (ny*nz)
          WRITE (18,'(e13.2)') obsvarx(i)
       END DO
       END IF

! ########## End i-level observation variance ##########

! ##### Calculates j-level observation variance #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          avgy(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                avgy(j) = avgy(j) + obs(i,j,k)
             END DO
          END DO
       END DO

       DO j=1,ny
          avgy(j) = avgy(j) / (nx * nz)
       END DO
          
       DO j=1,ny
          obsvary(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                obsvary(j) = obsvary(j) + (obs(i,j,k) - avgy(j))**2
             END DO
          END DO
       END DO
       
       WRITE (18, FMT = 802)
       
       DO j=1,ny
          obsvary(j) = obsvary(j) / (nx * nz)
          WRITE (18,'(e13.2)') obsvary(j)
       END DO
       END IF

! ########## End j-level observation variance ##########

! ##### Calculates k-level observation variance #####

       IF (klevel == 1) THEN
       DO k=1,nz
          avgz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                avgz(k) = avgz(k) + obs(i,j,k)
             END DO
          END DO
       END DO
       
       DO k=1,nz
          avgz(k) = avgz(k) / (nx * ny)
       END DO
          
       DO k=1,nz
          obsvarz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                obsvarz(k) = obsvarz(k) + (obs(i,j,k) - avgz(k))**2
             END DO
          END DO
       END DO
       
       WRITE (18, FMT = 803)
       
       DO k=1,nz
          obsvarz(k) = obsvarz(k) / (nx * ny)
          WRITE (18,'(e13.2)') obsvarz(k)
       END DO
       END IF

! ########## End k-level forecast variance ##########

       CLOSE (18)

800    FORMAT ("Observation Variance with all data"/)
801    FORMAT (/"Observation Variance with i-level dependence"/)
802    FORMAT (/"Observation Variance with j-level dependence"/)
803    FORMAT (/"Observation Variance with k-level dependence"/)
       CLOSE (UNIT = 18)


       END SUBROUTINE OBSERVATION_VARIANCE



!#######################################################################

       SUBROUTINE FORECAST_STD_DEVIATION (nx, ny, nz, model, ilevel, &
         jlevel, klevel)


!  Calculates the standard deviation for the hypothetical forecast data
!  set and writes this to the "FcstStdDev.txt" file
!
!  Formula:
!
!  Std Dev = SQRT[ SUM[ (x -x(bar))**2 / N]]
!
!  x is the value of each individual forecast
!  x(bar) is the mean of the set of forecasts
!  N is the number of forecasts

!#######################################################################

       IMPLICIT none
             
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num

       REAL :: stddev_all   
       REAL :: stddevz(nz)
       REAL :: stddevy(ny)
       REAL :: stddevx(nx)
       REAL :: avg_all
       REAL :: avgz(nz)
       REAL :: avgy(ny)
       REAL :: avgx(nx)
       REAL :: model(nx,ny,nz)

! ##### Calculates the standard deviation for all data #####

       avg_all = 0.0
       
       DO k=1,nz
          DO j=1,ny
             DO i=1,nx   
                avg_all = avg_all + model(i,j,k)
             END DO
          END DO
       END DO
       
       avg_all = avg_all / (nx * ny * nz)

       DO k=1,nz
          avgz(k) = avgz(k) / (nx * ny)
       END DO
       
       stddev_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                stddev_all = stddev_all + (model(i,j,k) - avg_all)**2
             END DO
          END DO
       END DO

       OPEN (UNIT = 20, FILE = "FcstStdDev.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
          
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF   
       
       WRITE (20, FMT = 100)
       
       stddev_all = stddev_all / (nx * ny * nz)
       stddev_all = SQRT(stddev_all)
       WRITE (20,'(e13.2)') stddev_all

! ########## End standard deviation for all data ##########

! ##### Calculates i-level standard deviation #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          avgx(i) = 0.0
          DO j=1,ny
             DO k=1,nz   
                avgx(i) = avgx(i) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO i=1,nx
          avgx(i) = avgx(i) / (ny * nz)
       END DO

       DO i=1,nx
          stddevx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                stddevx(i) = stddevx(i) + (model(i,j,k) - avgx(i))**2
             END DO
          END DO
       END DO
       
       WRITE (20, FMT = 101)
       
       DO i=1,nx
          stddevx(i) = stddevx(i) / (ny * nz)
          stddevx(i) = SQRT(stddevx(i))         
          WRITE (20,'(e13.2)') stddevx(i)
       END DO
       END IF

! ########## End i-level standard deviation ##########


! ##### Calculates j-level standard deviation #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          avgy(j) = 0.0
          DO k=1,nz
             DO i=1,nx   
                avgy(j) = avgy(j) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO j=1,ny
          avgy(j) = avgy(j) / (nx * nz)
       END DO
       
       DO j=1,ny
          stddevy(j) = 0.0
          DO k=1,nz
             DO i=1,nx
                stddevy(j) = stddevy(j) + (model(i,j,k) - avgy(j))**2
             END DO
          END DO
       END DO
       
       WRITE (20, FMT = 102)
       
       DO j=1,ny
          stddevy(j) = stddevy(j) / (nx * nz)
          stddevy(j) = SQRT(stddevy(j))
          WRITE (20,'(e13.2)') stddevy(j)
       END DO
       END IF

! ########## End j-level standard deviation ##########

! ##### Calculates k-level standard deviation #####


       IF (klevel == 1) THEN
       DO k=1,nz
          avgz(k) = 0.0
          DO j=1,ny
             DO i=1,nx   
                avgz(k) = avgz(k) + model(i,j,k)
             END DO
          END DO
       END DO
       
       DO k=1,nz
          avgz(k) = avgz(k) / (nx * ny)
       END DO
       
       DO k=1,nz
          stddevz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                stddevz(k) = stddevz(k) + (model(i,j,k) - avgz(k))**2
             END DO
          END DO
       END DO
       
       WRITE (20, FMT = 103)
       
       DO k=1,nz
          stddevz(k) = stddevz(k) / (nx * ny)
          stddevz(k) = SQRT(stddevz(k))
          WRITE (20,'(e13.2)') stddevz(k)
       END DO
       END IF

! ########## End k-level standard deviation ##########

       CLOSE(20)

100    FORMAT ("Forecast standard deviation for all data"/)
101    FORMAT (/"Forecast std dev with i-level dependence"/)
102    FORMAT (/"Forecast std dev with j-level dependence"/)
103    FORMAT (/"Forecast std dev with k-level dependence"/)


       END SUBROUTINE FORECAST_STD_DEVIATION





!#######################################################################

       SUBROUTINE OBSERVATION_STD_DEVIATION (nx, ny, nz, obs, ilevel, &
         jlevel, klevel)


!  Calculates the standard deviation for the hypothetical observation
!  set and writes this to the "ObsStdDev.txt" file
!
!  Formula:
!
!  Std Dev = SQRT[ SUM[ (x -x(bar))**2 / N]]
!
!  x is the value of each individual observation 
!  x(bar) is the mean of the set of observations 
!  N is the number of observations 

!#######################################################################

       IMPLICIT none
             
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num

       REAL :: stddev_all   
       REAL :: stddevz(nz)
       REAL :: stddevy(ny)
       REAL :: stddevx(nx)
       REAL :: avg_all
       REAL :: avgz(nz)
       REAL :: avgy(ny)
       REAL :: avgx(nx)
       REAL :: obs(nx,ny,nz)

! ##### Calculates the standard deviation for all data #####

       avg_all = 0.0
       
       DO k=1,nz
          DO j=1,ny
             DO i=1,nx   
                avg_all = avg_all + obs(i,j,k)
             END DO
          END DO
       END DO
       
       avg_all = avg_all / (nx * ny * nz)

       stddev_all = 0.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                stddev_all = stddev_all + (obs(i,j,k) - avg_all)**2
             END DO
          END DO
       END DO

       OPEN (UNIT = 21, FILE = "ObsStdDev.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
       
       WRITE (21, FMT = 200)
       
       stddev_all = stddev_all / (nx * ny * nz)
       stddev_all = SQRT(stddev_all)
       WRITE (21,'(e13.2)') stddev_all

! ########## End standard deviation for all data ##########

! ##### Calculates i-level standard deviation #####

       IF (ilevel == 1) THEN
       DO i=1,nx
          avgx(i) = 0.0
          DO j=1,ny
             DO k=1,nz   
                avgx(i) = avgx(i) + obs(i,j,k)
             END DO
          END DO
       END DO
       
       DO i=1,nx
          avgx(i) = avgx(i) / (ny * nz)
       END DO

       DO i=1,nx
          stddevx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                stddevx(i) = stddevx(i) + (obs(i,j,k) - avgx(i))**2
             END DO
          END DO
       END DO
       
       WRITE (21, FMT = 201)
       
       DO i=1,nx
          stddevx(i) = stddevx(i) / (ny * nz)
          stddevx(i) = SQRT(stddevx(i))         
          WRITE (21,'(e13.2)') stddevx(i)
       END DO
       END IF

! ########## End i-level standard deviation ##########


! ##### Calculates j-level standard deviation #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          avgy(j) = 0.0
          DO k=1,nz
             DO i=1,nx   
                avgy(j) = avgy(j) + obs(i,j,k)
             END DO
          END DO
       END DO
       
       DO j=1,ny
          avgy(j) = avgy(j) / (nx * nz)
       END DO
       
       DO j=1,ny
          stddevy(j) = 0.0
          DO k=1,nz
             DO i=1,nx
                stddevy(j) = stddevy(j) + (obs(i,j,k) - avgy(j))**2
             END DO
          END DO
       END DO
       
       WRITE (21, FMT = 202)
       
       DO j=1,ny
          stddevy(j) = stddevy(j) / (nx * nz)
          stddevy(j) = SQRT(stddevy(j))
          WRITE (21,'(e13.2)') stddevy(j)
       END DO
       END IF

! ########## End j-level standard deviation ##########

! ##### Calculates k-level standard deviation #####


       IF (klevel == 1) THEN
       DO k=1,nz
          avgz(k) = 0.0
          DO j=1,ny
             DO i=1,nx   
                avgz(k) = avgz(k) + obs(i,j,k)
             END DO
          END DO
       END DO
       
       DO k=1,nz
          avgz(k) = avgz(k) / (nx * ny)
       END DO
       
       DO k=1,nz
          stddevz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                stddevz(k) = stddevz(k) + (obs(i,j,k) - avgz(k))**2
             END DO
          END DO
       END DO
       
       WRITE (21, FMT = 203)
       
       DO k=1,nz
          stddevz(k) = stddevz(k) / (nx * ny)
          stddevz(k) = SQRT(stddevz(k))
          WRITE (21,'(e13.2)') stddevz(k)
       END DO
       END IF

! ########## End k-level standard deviation ##########

       CLOSE(21)

200    FORMAT ("Observation standard deviation for all data"/)
201    FORMAT (/"Observation std dev with i-level dependence"/)
202    FORMAT (/"Observation std dev with j-level dependence"/)
203    FORMAT (/"Observation std dev with k-level dependence"/)


       END SUBROUTINE OBSERVATION_STD_DEVIATION




!#######################################################################

       SUBROUTINE BIAS_CALC (nx, ny, nz, data, model, ilevel, jlevel, &
         klevel,bias_all)

       
!  Calculates the bias along the k level of the grid between the
!  hypothetical data set and model forecast and writes this
!  to the "JBias.txt" file
!
!  Formula:
!  
!  Bias = (1/n) * SUM[k=1 to n] (y(k) - o(k))
!      
!  n is the number of forecast / observation pairs
!  y(k) is the kth of n forecasts
!  o(k) is the kth of n observations

!####################################################################### 


       IMPLICIT none
       
       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios, num
       INTEGER :: biascount
       
       REAL :: bias_all
       REAL :: biasz(nz)
       REAL :: biasy(ny)
       REAL :: biasx(nx)
       REAL :: data(nx,ny,nz), model(nx,ny,nz)
       
! ##### Calclulates the bias for all data #####

       bias_all = 0.0
       biascount=0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                IF (data(i,j,k).ne.-99.9) THEN
                bias_all = bias_all + (model(i,j,k)-data(i,j,k))
                biascount=biascount+1
                END IF
             END DO
          END DO
       END DO

       IF ( biascount == 0) THEN
          bias_all = -99.9
          RETURN
       END IF

!       OPEN (UNIT = 16, FILE = "Bias.txt", IOSTAT = ios, &
!       FORM = "FORMATTED")

!       IF (ios /= 0) THEN
!          PRINT*, "Error opening file"
!       END IF

!       WRITE (16, FMT = 900)

!       bias_all = bias_all / (nx * ny * nz)
        bias_all = bias_all / (biascount)
!       WRITE (16,'(f7.2)') bias_all

! ########## End bias for all data ##########

! ##### Calculates i-level bias ######

       IF (ilevel == 1) THEN
       DO i=1,nx
          biasx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                biasx(i) = biasx(i) + (model(i,j,k)-data(i,j,k))
             END DO
          END DO
       END DO
       
       WRITE (16, FMT = 901)
       
       DO i = 1,nx
          biasx(i) = biasx(i) / (ny * nz)
          WRITE (16,'(f7.2)') biasx(i)
       END DO
       END IF

! ########## End i-level bias ##########

! ##### Calculates j-level bias #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          biasy(j) = 0.0
          DO k=1,nz
             DO i=1,nx
                biasy(j) = biasy(j) + (model(i,j,k)-data(i,j,k))
             END DO
          END DO
       END DO
       
       WRITE (16, FMT = 902)
       
       DO j = 1,ny
          biasy(j) = biasy(j) / (nx * nz)
          WRITE (16,'(f7.2)') biasy(j)
       END DO
       END IF

! ########## End j-level bias ##########

! ##### Calculate k-level bias #####

       IF (klevel == 1) THEN
       DO k=1,nz
          biasz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                biasz(k) = biasz(k) + (model(i,j,k)-data(i,j,k))
             END DO
          END DO
       END DO
       
       WRITE (16, FMT = 903)

       DO k = 1,nz
          biasz(k) = biasz(k) / (nx * ny)
          WRITE (16,'(f7.2)') biasz(k)
       END DO
       END IF

! ########## End k-level bias ##########

       
       CLOSE (UNIT = 16)

900    FORMAT ("Bias for all data"/)
901    FORMAT (/"Bias with i-level dependence"/)
902    FORMAT (/"Bias with j-level dependence"/)
903    FORMAT (/"Bias with k-level dependence"/)


       END SUBROUTINE BIAS_CALC




!#######################################################################
       
       SUBROUTINE THREAT_SCORE (noe, fcstprob, obsevent)
       
       
!  Calculates the threat score using hypothetical data from a 
!  probability forecast and writes this calculation to the
!  "TS.txt" file.  This calculation also goes by the name
!  Critical Success Index (CSI).
!
!  Formula:
!
!  Threat Score =  CFE(y) / [ FE(y) + OE(y) - CFE(y) ]
!      
!  
!  CFE(y) is the number of correctly forecast 'yes' events
!  FE(y) is the number of forecast 'yes' events
!  OE(y) is the number of observed 'yes' events

!#######################################################################

       IMPLICIT none
       
       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
       REAL :: cfey  !  number of correctly forecast 'yes' events
       REAL :: fey   !  number of forecast 'yes' events
       REAL :: oey   !  number of observed 'yes' events
          
       REAL :: threat   !  threat score
       
       REAL :: fcstprob(noe), obsevent(noe)

       cfey = 0
       fey  = 0
       oey  = 0

       DO g=1,noe
          IF ((fcstprob(g).gt.(.50)).AND.(obsevent(g) == 1.0)) THEN
             cfey = cfey + 1
          END IF
       END DO

       DO g=1,noe
          IF (fcstprob(g).gt.(.50)) THEN
             fey = fey + 1
          END IF
       END DO

       DO g=1,noe
          IF (obsevent(g) == 1) THEN
             oey = oey + 1
          END IF
       END DO

       threat = cfey / ( fey + oey - cfey )

       OPEN (UNIT = 23, FILE = "TS.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN   
          PRINT*, "Error opening file"
       END IF
       
       WRITE (23, FMT = 100)
       WRITE (23, '(f7.4)') threat

       CLOSE (23)

100    FORMAT ("Threat Score"/)


       END SUBROUTINE THREAT_SCORE



!#######################################################################
       
       SUBROUTINE EQUITABLE_THREAT_SCORE (noe, fcstprob, obsevent)
       
       
!  Calculates the equitable threat score using hypothetical data from a
!  probability forecast and writes this calculation to the
!  "ETS.txt" file. 
!      
!  Formula:
!
!  Equitable Threat Score =(CFE(y) - RH) / [FE(y) + OE(y) - CFE(y) - RH]
!
!  RH = (FE(y) * OE(y)) / n
!
!  RH is the number of hits above the threshold that would be expected
!     in a random forecast
!  CFE(y) is the number of correctly forecast 'yes' events
!  FE(y) is the number of forecast 'yes' events
!  OE(y) is the number of observed 'yes' events
!  n is the total number of forecasts
       
!#######################################################################
       
       IMPLICIT none
             
       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
       REAL :: cfey  !  number of correctly forecast 'yes' events
       REAL :: fey   !  number of forecast 'yes' events
       REAL :: oey   !  number of observed 'yes' events
       REAL :: rh    !  random hits             

       REAL :: ets   !  equitable threat score
       
       REAL :: fcstprob(noe), obsevent(noe)
       
       cfey = 0
       fey  = 0
       oey  = 0

       DO g=1,noe
          IF ((fcstprob(g).gt.(.50)).AND.(obsevent(g) == 1.0)) THEN
             cfey = cfey + 1
          END IF
       END DO
       
       DO g=1,noe
          IF (fcstprob(g).gt.(.50)) THEN
             fey = fey + 1
          END IF
       END DO

       DO g=1,noe
          IF (obsevent(g) == 1) THEN
             oey = oey + 1
          END IF
       END DO

       rh = (fey * oey) / noe


       ets = (cfey - rh) / (fey + oey - cfey -rh)
       
       OPEN (UNIT = 24, FILE = "ETS.txt", IOSTAT = ios, &
       FORM = "FORMATTED")

       IF (ios /= 0) THEN 
          PRINT*, "Error opening file"
       END IF
 
       WRITE (24, FMT = 100)
       WRITE (24, '(f7.4)') ets

       CLOSE (24)

100    FORMAT ("Equitable Threat Score"/)


       END SUBROUTINE EQUITABLE_THREAT_SCORE




!#######################################################################

       SUBROUTINE PROB_OF_DETECTION (noe, fcstprob, obsevent)
       
       
!  Calculates the probability of detection using hypothetical data 
!  from a probability forecast and writes this calculation to the
!  "POD.txt" file.
!
!  Formula:
!
!  Probability of Detection = CFE(y) / OE(y)
!
!  CFE(y) is the number of correctly forecast 'yes' events
!  OE(y) is the number of observed 'yes' events
          
!#######################################################################
          
       IMPLICIT none

       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
       REAL :: cfey  !  number of correctly forecast 'yes' events
       REAL :: oey   !  number of observed 'yes' events
       
       REAL :: pod   !  probability of detection
             
       REAL :: fcstprob(noe), obsevent(noe)
       
       cfey = 0
       oey  = 0

       DO g=1,noe
          IF ((fcstprob(g).gt.(.50)).AND.(obsevent(g) == 1.0)) THEN
             cfey = cfey + 1
          END IF
       END DO
       
       DO g=1,noe
          IF (obsevent(g) == 1) THEN
             oey = oey + 1
          END IF
       END DO

       pod = cfey / oey

       OPEN (UNIT = 25, FILE = "POD.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF
             
       WRITE (25, FMT = 100)
       WRITE (25, '(f7.4)') pod
   
       CLOSE (25)
          
100    FORMAT ("Probability of Detection"/)
          
       
       END SUBROUTINE PROB_OF_DETECTION




!#######################################################################
       
       SUBROUTINE FALSE_ALARM_RATE (noe, fcstprob, obsevent)
       
       
!  Calculates the false alarm rate using hypothetical data
!  from a probability forecast and writes this calculation to the
!  "FAR.txt" file.
!
!  Formula:
!            
!  False Alarm Rate = IFE(n) / FE(y)
!
!  IFE(n) is the number of incorrectly forecast 'no' events
!     (ie: The event was forecast but not observed.)
!  FE(y) is the number of forecast 'yes' events
       
!#######################################################################
       
       IMPLICIT none
             
       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
       REAL :: ifen  !  number of incorrectly forecast 'yes' events
       REAL :: fey   !  number of forecast 'yes' events

       REAL :: far   !  false alarm rate
     
       REAL :: fcstprob(noe), obsevent(noe)
       
       ifen = 0
       fey  = 0

       DO g=1,noe
          IF ((fcstprob(g).gt.(.50)).AND.(obsevent(g) == 0.0)) THEN
             ifen = ifen + 1
          END IF 
       END DO
       
       DO g=1,noe
          IF (fcstprob(g).gt.(.50)) THEN
             fey = fey + 1
          END IF
       END DO

       far = ifen / fey

       OPEN (UNIT = 26, FILE = "FAR.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF  

       WRITE (26, FMT = 100)
       WRITE (26, '(f7.4)') far
             
       CLOSE (26)
       
100    FORMAT ("False Alarm Rate"/)
       
          
       END SUBROUTINE FALSE_ALARM_RATE





!######################################################################

       SUBROUTINE HIT_RATE (noe, fcstprob, obsevent)
       

!  Calculates the hit rate using hypothetical data from a probability
!  forecast and writes this calculation to the "HR.txt" file.
!      
!  Formula:
!
!  Hit Rate = (CFE(y) + CFE(n)) / n
!
!  CFE(y) is the number of correctly forecast 'yes' events
!  CFE(n) is the number of correctly forecast 'no' events
!  n is the total number of forecasts
       
!#######################################################################
       
       IMPLICIT none
             
       INTEGER :: noe   !  number of events
       INTEGER :: g, ios
       REAL :: cfey  !  number of correctly forecast 'yes' events
       REAL :: cfen  !  number of correctly forecast 'no' events
       
       REAL :: hr   !  hit rate
     
       REAL :: fcstprob(noe), obsevent(noe)
       
       cfey = 0.0
       cfen = 0.0
       
       DO g=1,noe
          IF ((fcstprob(g).gt.(.50)).AND.(obsevent(g) == 1.0)) THEN
             cfey = cfey + 1
          END IF 
       END DO

       DO g=1,noe
          IF ((fcstprob(g).le.(.50)).AND.(obsevent(g) == 0.0)) THEN
             cfen = cfen + 1
          END IF
       END DO

       hr = (cfey + cfen) / noe

       OPEN (UNIT = 27, FILE = "HR.txt", IOSTAT = ios, &
       FORM = "FORMATTED")
       
       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF

       WRITE (27, FMT = 100)
       WRITE (27, '(f7.4)') hr
 
       CLOSE (27)
 
100    FORMAT ("Hit Rate"/)


       END SUBROUTINE HIT_RATE



!#######################################################################

       SUBROUTINE SKILL_SCORE (nx, ny, nz, model, ilevel, &
         jlevel, klevel)


!  Calculates the skill score using the hypothetical model forecast
!  and a hypothetical reference value writes this to the "SS.txt" file
!
!  The file "SS.txt" will always contain the mean absolute error for 
!  the entire data set at the top of the file.  It will also contain
!  the mean absolute error based on the i, j, or k levels depending
!  upon which of those options are chosen in the "verif.input" file.
!
!
!  Formula:
!
!  MAE = (1/n) * SUM (forecast - reference) / (1 - reference)
!
!  n is the number of forecasts
!  forecast is the forecast value
!  reference is the "climotology" or "persistance" value               

!#######################################################################

       IMPLICIT none

       INTEGER :: nx   !  number of x grid points
       INTEGER :: ny   !  number of y grid points
       INTEGER :: nz   !  number of z grid points
       INTEGER :: ilevel   !  toggles ilevel calculation
       INTEGER :: jlevel   !  toggles jlevel calculation
       INTEGER :: klevel   !  toggles klevel calculation
       INTEGER :: i, j ,k, ios

       REAL :: ssz(nz)
       REAL :: ssy(ny)
       REAL :: ssx(nx)
       REAL :: ss_all   !  variable for SS of entire data set 
       REAL :: num
       REAL :: ref   !  reference value
       REAL :: model(nx,ny,nz)

!  ##### Calculates skill score for entire data set mae_all #####

       ss_all = 0.0
       ref = 50.0

       DO k=1,nz
          DO j=1,ny
             DO i=1,nx
                num = 0.0
                num = (model(i,j,k) - ref) / (1 - ref) 
                ss_all = ss_all + num
             END DO
          END DO
       END DO

       OPEN (UNIT = 28, FILE = "SS.txt", IOSTAT = ios,  &
       FORM = "FORMATTED")

       IF (ios /= 0) THEN
          PRINT*, "Error opening file"
       END IF

       WRITE (28, FMT = 200)

       ss_all = ss_all / (nx * ny * nz)
       WRITE (28,'(f7.2)') ss_all

! ########## End of ss_all calculation ##########


! ##### Calculates i-level skill score, ssx ##### 

       IF (ilevel == 1) THEN
       DO i=1,nx
          ssx(i) = 0.0
          DO j=1,ny
             DO k=1,nz
                num = 0.0
                num = (model(i,j,k) - ref) / (1 - ref)
                ssx(i) = ssx(i) + num
             END DO
          END DO
       END DO
       
       WRITE (28, FMT = 201)

       DO i=1,nx   
          ssx(i) = ssx(i) / (ny * nz)
          WRITE (28,'(f7.2)') ssx(i)
       END DO
       END IF

! ########## End of i-level skill score ##########

! ##### Calculates j-level skill score, ssy #####

       IF (jlevel == 1) THEN
       DO j=1,ny
          ssy(j) = 0.0
          DO i=1,nx
             DO k=1,nz
                num = 0.0
                num = (model(i,j,k) - ref) / (1 - ref)
                ssy(j) = ssy(j) + num
             END DO
          END DO
       END DO

       WRITE (28, FMT = 202)
       
       DO j=1,ny 
          ssy(j) = ssy(j) / (nx * nz)
          WRITE (28,'(f7.2)') ssy(j)
       END DO
       END IF

! ######### End of j-level skill score ##########

! ##### Calculates the k-level skill score #####

       IF (klevel == 1) THEN
       DO k=1,nz
          ssz(k) = 0.0
          DO j=1,ny
             DO i=1,nx
                num = 0.0
                num = (model(i,j,k) - ref) / (1 - ref)
                ssz(k) = ssz(k) + num
             END DO
          END DO
       END DO

       WRITE (28, FMT = 203)
       
       DO k=1,nz 
          ssz(k) = ssz(k) / (nx * ny)
          WRITE (28,'(f7.2)') ssz(k)
       END DO
       END IF

! ######### End of k-level skill score #########

       CLOSE (UNIT = 28)


200    FORMAT ("Skill Score for all data"/)
201    FORMAT (/"Skill Score with i-level dependence"/)
202    FORMAT (/"Skill Score with j-level dependence"/)
203    FORMAT (/"Skill Score with k-level dependence"/)

       END SUBROUTINE SKILL_SCORE
