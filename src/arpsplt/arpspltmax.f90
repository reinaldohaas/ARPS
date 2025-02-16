PROGRAM pltmax
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                PROGRAM PLTMAX                        ######
!######                                                      ######
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
!  This program produces graphic plots of time series of
!  the domain maximum and minimum of u, v, w, potential temperature
!  perturbation, pressure perturbation and the water quantities.
!
!  A max./min. data file produced by the ARPS model run is
!  used as the input data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue, CAPS/OU.
!  6/26/1992
!
!  MODIFICATION HISTORY:
!
!  4/13/93 (K. Droegemeier)
!  Added reading and plotting part for reflectivity, vorticity etc.
!
!  8/13/95 (A. Shapiro)
!  Minor documentation clean-up.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!
!  CALCULATED DATA ARRAYS:
!
!
!  WORKING ARRAYS:
!
!    tem1     Temporary working array.
!    tem2     Temporary working array.
!    tem3     Temporary working array.
!    tem4     Temporary working array.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nmax
  PARAMETER (nmax = 1000)  ! Maxmimum array size
!
  REAL :: time (nmax)
!
  REAL :: umax (nmax)
  REAL :: umin (nmax)
  REAL :: vmax (nmax)
  REAL :: vmin (nmax)
  REAL :: wmax (nmax)
  REAL :: wmin (nmax)
  REAL :: ptmax(nmax)
  REAL :: ptmin(nmax)
  REAL :: pmax (nmax)
  REAL :: pmin (nmax)

  REAL :: qvmax(nmax)
  REAL :: qcmax(nmax)
  REAL :: qrmax(nmax)
  REAL :: qimax(nmax)
  REAL :: qsmax(nmax)
  REAL :: qhmax(nmax)

  REAL :: ulmax(nmax)
  REAL :: vilmax(nmax)
  REAL :: vorlmax(nmax)
  REAL :: vorumax(nmax)
  REAL :: sfcrain(nmax)
  REAL :: refmax(nmax)
  REAL :: zrefmax(nmax)

  REAL :: uulim,ullim   ! upper and lower limit of y-axis for u graph
  REAL :: vulim,vllim   ! upper and lower limit of y-axis for v graph
  REAL :: wulim,wllim   ! upper and lower limit of y-axis for w graph
  REAL :: ptllim,ptulim ! upper and lower limit of y-axis for pt graph
  REAL :: pulim,pllim   ! upper and lower limit of y-axis for p graph
  REAL :: qulim,qllim   ! upper and lower limit of y-axis for q graph

  REAL :: ululim,ulllim ! upper and lower limit of y-axis for max wind graph
  REAL :: vorlulim,vorlllim ! upper and lower limits for vort below
  REAL :: voruulim,vorullim ! upper and lower limits for vort above
  REAL :: refulim,refllim   ! upper and lower limits for reflectivity

  CHARACTER (LEN=80)  :: runname
  CHARACTER (LEN=256) :: filename
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, npoints, nxpic,nypic,nchanl,lenstr
  REAL :: x1,x2,y1,y2,xlimit,ylimit,angl,timmin,timmax

!-----------------------------------------------------------------------
!
!  Plotting control parameters entered by user
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256) :: outfilename
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
!  Obtain user-specified control parameters for reading in
!  max/min data and plotting time series.
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(/ 16(/5x,a)//)')                                            &
     '###############################################################', &
     '###############################################################', &
     '#####                                                     #####', &
     '#####               Welcome to ARPSPLTMAX                 #####', &
     '#####                                                     #####', &
     '#####           A program that plots the graphs           #####', &
     '#####            of the max./min. time series             #####', &
     '#####                                                     #####', &
     '#####           The graphic plotting is based             #####', &
     '#####             on graphic package ZXPLOT               #####', &
     '#####                                                     #####', &
     '#####              By Ming Xue CAPS/OU                    #####', &
     '#####                                                     #####', &
     '###############################################################', &
     '###############################################################'

!
!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!

  filename = ' '
  WRITE(6,'(/a/)')'Please give the name of the input data file:'
  READ(5,*) filename

  lenstr = LEN_TRIM(filename)

  WRITE(6,'(/a,a)')' The data set name was ', filename(1:lenstr)

  nchanl = 63

  OPEN(UNIT=nchanl,FILE=filename(1:lenstr),STATUS='old',FORM='formatted')

  WRITE(6,'(/a,a,a,i3/)')' File ', filename(1:lenstr),                  &
       ' was opened using Fortran unit ',nchanl


!
!-----------------------------------------------------------------------
!
!  Read in the upper and lower limit of the y axis for each
!  variable.
!
!-----------------------------------------------------------------------
!

  WRITE(6,'(/a/)')                                                      &
      'Please specify the lower limit of y axis of the umax/umin plot.'
  READ(5,*) ullim
  WRITE(6,*) 'Input was ',ullim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the upper limit of y axis of the umax/umin plot.'
  READ(5,*) uulim
  WRITE(6,*) 'Input was ',uulim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the lower limit of y axis of the vmax/vmin plot.'
  READ(5,*) vllim
  WRITE(6,*) 'Input was ',vllim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the upper limit of y axis of the vmax/vmin plot.'
  READ(5,*) vulim
  WRITE(6,*) 'Input was ',vulim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the lower limit of y axis of the wmax/wmin plot.'
  READ(5,*) wllim
  WRITE(6,*) 'Input was ',wllim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the upper limit of y axis of the wmax/wmin plot.'
  READ(5,*) wulim
  WRITE(6,*) 'Input was ',wulim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the lower limit of y axis of ',                   &
      'the ptmax/ptmin plot.'
  READ(5,*) ptllim
  WRITE(6,*) 'Input was ',ptllim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the upper limit of y axis of ',                   &
      'the ptmax/ptmin plot.'
  READ(5,*) ptulim
  WRITE(6,*) 'Input was ',ptulim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the lower limit of y axis of the pmax/pmin plot.'
  READ(5,*) pllim
  WRITE(6,*) 'Input was ',pllim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the upper limit of y axis of the pmax/pmin plot.'
  READ(5,*) pulim
  WRITE(6,*) 'Input was ',pulim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the lower limit of y axis of the qmax/qmin plot.'
  READ(5,*) qllim
  WRITE(6,*) 'Input was ',qllim

  WRITE(6,'(/a/)')                                                      &
      'Please specify the upper limit of y axis of the qmax/qmin plot.'
  READ(5,*) qulim
  WRITE(6,*) 'Input was ',qulim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the lower limit of y axis of ',                   &
      'the low-level wind speed plot'
  READ(5,*) ulllim
  WRITE(6,*) 'Input was ',ulllim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the upper limit of y axis of ',                   &
      'the low-level wind speed plot'
  READ(5,*) ululim
  WRITE(6,*) 'Input was ',ululim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the lower limit of y axis of ',                   &
      'the vertical vorticity below specified height'
  READ(5,*) vorlllim
  WRITE(6,*) 'Input was ',vorlllim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the upper limit of y axis of ',                   &
      'the vertical vorticity below specified height'
  READ(5,*) vorlulim
  WRITE(6,*) 'Input was ',vorlulim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the lower limit of y axis of ',                   &
      'the vertical vorticity above specified height'
  READ(5,*) vorullim
  WRITE(6,*) 'Input was ',vorullim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the upper limit of y axis of ',                   &
      'the vertical vorticity above specified height'
  READ(5,*) voruulim
  WRITE(6,*) 'Input was ',voruulim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the lower limit of y axis of ',                   &
      'the max radar reflectivity'
  READ(5,*) refllim
  WRITE(6,*) 'Input was ',refllim

  WRITE(6,'(/a,a/)')                                                    &
      'Please specify the upper limit of y axis of ',                   &
      'the max radar reflectivity'
  READ(5,*) refulim
  WRITE(6,*) 'Input was ',refulim

  WRITE(6,'(/a/)')                                                      &
      'How many columns of plots do you want in each frame?'
  READ(5,*) nxpic
  WRITE(6,*) 'Input was ',nxpic

  WRITE(6,'(/a/)')                                                      &
      'How many raws of plots do you want in each frame?'
  READ(5,*) nypic
  WRITE(6,*) 'Input was ',nypic

!
!-----------------------------------------------------------------------
!
!  Read in the max./min. data
!
!-----------------------------------------------------------------------
!

  READ(nchanl,*,ERR=999,END=992) runname

  WRITE(6,'(/a/)') 'Reading data from ',runname

  READ(nchanl,*)

  DO i = 1,nmax

    READ(nchanl,*,ERR=999,END=992)                                      &
        time(i),umin(i),umax(i),vmin(i),vmax(i),wmin(i),wmax(i),        &
        ptmin(i),ptmax(i),pmin(i),pmax(i),qvmax(i),qcmax(i),qrmax(i),   &
        qimax(i),qsmax(i),qhmax(i),                                     &
        ulmax(i),vorlmax(i),vorumax(i),refmax(i)

    pmin(i)=pmin(i)*0.01
    pmax(i)=pmax(i)*0.01
    vorlmax(i)=vorlmax(i)*1000.
    vorumax(i)=vorumax(i)*1000.
    time(i) = time(i)/60.0

  END DO

  992   CONTINUE

  npoints = i-1   ! The number of data samples in each time series

  timmin = time(1)
  timmax = time(npoints)
!
!-----------------------------------------------------------------------
!
!  Initialize plotting package and set plotting space parameters
!
!-----------------------------------------------------------------------
!
  CALL xdevic_new(1,outfilename,0,0)

  IF( nxpic == 1 .AND. nypic == 1) THEN
    CALL xdspac(0.80)
  ELSE
    CALL xdspac(0.90)
  END IF

  angl = 00.0
!
!-----------------------------------------------------------------------
!
!  Set page layout parameters
!
!-----------------------------------------------------------------------
!
  CALL xspace(nxpic, nypic, angl , xlimit,ylimit)
!
!-----------------------------------------------------------------------
!
!  Set the margins between pictures in nondimensional units.
!
!-----------------------------------------------------------------------
!
  CALL xpmagn(0.05,0.05)

!
!-----------------------------------------------------------------------
!
!  Plot umax, umin curves
!
!-----------------------------------------------------------------------
!
  CALL xaxfmt( '(i4)' )

  CALL xnwpic

  CALL xmap(timmin,timmax,ullim, uulim)
  CALL xbordr

  PRINT*, timmin, timmax

!  call xqmap(x1,x2,y1,y2)
!  call xaxnsz( 0.025*(y2-y1) )

  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Umin, Umax (m/s)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,umax,npoints)
  CALL xdash
  CALL xgraph(time,umin,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max U Plot Finished'
!
!-----------------------------------------------------------------------
!
!  Plot vmax, vmin curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,vllim, vulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Vmin, Vmax (m/s)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,vmax,npoints)
  CALL xdash
  CALL xgraph(time,vmin,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max V Plot Finished'

!
!-----------------------------------------------------------------------
!
!  Plot wmax, wmin curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,wllim,wulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Wmin, Wmax (m/s)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,wmax,npoints)
  CALL xdash
  CALL xgraph(time,wmin,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max W Plot Finished'

!
!-----------------------------------------------------------------------
!
!  Plot ptmax, ptmin curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,ptllim,ptulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'PTmin, PTmax (K)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,ptmax,npoints)
  CALL xdash
  CALL xgraph(time,ptmin,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max PT Plot Finished'
!
!-----------------------------------------------------------------------
!
!  Plot pmax, pmin curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,pllim, pulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Pmin, Pmax (mb)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,pmax,npoints)
  CALL xdash
  CALL xgraph(time,pmin,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max P Plot Finished'

!
!-----------------------------------------------------------------------
!
!  Plot qvmax, qcmax and qrmax curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,qllim,qulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5,                            &
             'Qvmax, Qcmax, Qrmax (kg/kg)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )
  CALL xfull
  CALL xgraph(time,qvmax,npoints)
  CALL xdash
  CALL xgraph(time,qcmax,npoints)
  CALL runlab( runname )
  CALL xbrokn(8,4,8,1)
  CALL xgraph(time,qrmax,npoints)
  CALL xfull

  PRINT *, 'Max Water Plot Finished'

!
!-----------------------------------------------------------------------
!
!  Plot qsmax, qimax and qhmax curves
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,qllim,qulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5,                            &
             'Qsmax, Qimax, Qhmax (kg/kg)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )
  CALL xfull
  CALL xgraph(time,qsmax,npoints)
  CALL xdash
  CALL xgraph(time,qimax,npoints)

  CALL runlab( runname )
  CALL xbrokn(8,4,8,1)
  CALL xgraph(time,qhmax,npoints)
  CALL xfull
  PRINT *, 'Max Ice Plot Finished'
!
!-----------------------------------------------------------------------
!
!  Plot max low-level wind speed
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,ulllim,ululim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Sfc Wind Max (m/s)')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,ulmax,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max Low-Level Wind Plot Finished'
!
!-----------------------------------------------------------------------
!
!  Plot max vertical vorticity below specified height
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,vorlllim,vorlulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc(x1-0.10*(x2-x1),(y1+y2)*0.5,'Max Zeta(10**3) < 2 km')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,vorlmax,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max Vert Vorticity Below Plot Finished'

!
!-----------------------------------------------------------------------
!
!  Plot max vertical vorticity above specified height
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,vorullim,voruulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc(x1-0.10*(x2-x1),(y1+y2)*0.5,'Max Zeta(10**3) > 2 km')
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,vorumax,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max Vert Vorticity Above Plot Finished'
!
!-----------------------------------------------------------------------
!
!  Plot max radar reflectivity
!
!-----------------------------------------------------------------------
!
  CALL xnwpic
  CALL xmap(timmin,timmax,refllim,refulim)
  CALL xbordr
  CALL xaxes(timmin,0.0,0.0,0.0)
  CALL xqmap(x1,x2,y1,y2)
  CALL xchsiz(0.03*(y2-y1))
  CALL xchori(90.0)
  CALL xcharc( x1-0.10*(x2-x1), (y1+y2)*0.5, 'Max dBz' )
  CALL xchori( 0.0)
  CALL xcharc( (x1+x2)*0.5, y1-(y2-y1)*0.10, 'Time (min.)' )

  CALL xfull
  CALL xgraph(time,refmax,npoints)
  CALL xfull

  CALL runlab( runname )
  PRINT *, 'Max Reflectivity Plot Finished'


  CALL xgrend

  STOP

  999   CONTINUE

  WRITE(6,'(2(1x,a))')'Error encountered when reading input data,',     &
        'Job aborted.'
  STOP

END PROGRAM pltmax
!

SUBROUTINE runlab_1(runname)
!
!-----------------------------------------------------------------------
!
!  To plot a run label at the lower left cornor of the picture frame.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER (LEN=*) :: runname
  REAL :: xl,xr,yb,yt,rotang,xlimit,ylimit
  INTEGER :: nopic,nxpic,nypic

  CALL xqmap(xl,xr,yb,yt)
  CALL xqnpic(nopic)
  CALL xqspac(nxpic, nypic, rotang, xlimit, ylimit)

  IF( rotang == 0.0 ) THEN

    IF(nopic == nxpic*nypic -(nxpic-1)) THEN
      CALL xcharl( xl, yb-0.15*(yt-yb), runname )
    END IF

  ELSE

    IF(nopic == nypic*nxpic -(nypic-1)) THEN
      CALL xcharl( xl, yb-0.15*(yt-yb), runname )
    END IF

  END IF

  RETURN
END SUBROUTINE runlab_1
