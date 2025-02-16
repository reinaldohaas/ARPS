!########################################################################
!########################################################################
!######                                                            ######
!######                    MODULE EXTDIMS2                         ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

MODULE extdims2

!-----------------------------------------------------------------------
!
! External grid variables
!
!-----------------------------------------------------------------------
!
! MODIFICATION HISTORY:
! Eric Kemp, November 1999
! Added unique variables for GRIB and GEMPAK files.  NOTE:  Currently
! this module is used only by CVT2VERIF; EXT2ARPS still uses the
! original EXTDIMS module.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
!
! RUC from NMC is 81x62x25
! RUC GRIB #87
!
!-----------------------------------------------------------------------

  INTEGER :: ruc87nx,ruc87ny,ruc87nz
  PARAMETER (ruc87nx=81,ruc87ny=62,ruc87nz=25)

!-----------------------------------------------------------------------
!
! ETA (GRIB #212, 40km) from NMC is 185x129x39
!
!-----------------------------------------------------------------------

  INTEGER :: eta212nx,eta212ny,eta212nz
  PARAMETER( eta212nx = 185, eta212ny = 129, eta212nz = 39 )

!-----------------------------------------------------------------------
!
! ETA (GEMPAK #104, 80km) from NMC is 147x110x39
!
!-----------------------------------------------------------------------

  INTEGER :: etagemnx,etagemny,etagemnz
  PARAMETER (etagemnx=147,etagemny=110,etagemnz=39)

!-----------------------------------------------------------------------
!
! RUC from Doplight is 81x62x41
! RUC GEMPAK
!
!-----------------------------------------------------------------------

  INTEGER :: rucgemnx,rucgemny,rucgemnz
  PARAMETER (rucgemnx=81,rucgemny=62,rucgemnz=41)

!-----------------------------------------------------------------------
!
! RUC-2 (GEMPAK) from Doplight is 151x113x41
!
!-----------------------------------------------------------------------

  INTEGER :: ruc2gemnx,ruc2gemny,ruc2gemnz
  PARAMETER (ruc2gemnx=151,ruc2gemny=113,ruc2gemnz=41)

!-----------------------------------------------------------------------
!
! RUC from Rossby is 93x65x19
!
!-----------------------------------------------------------------------

  INTEGER :: ruc211nx,ruc211ny,ruc211nz
  PARAMETER (ruc211nx=93,ruc211ny=65,ruc211nz=19)

!-----------------------------------------------------------------------
!
! OLAPS for 95 is 91x73x41
! NOTE:  OLAPS dimensions are passed to CVT2VERIF by the input file.
!
!-----------------------------------------------------------------------
!
! parameter (nx_ext=91,ny_ext=73,nz_ext=41)
!  
!-----------------------------------------------------------------------
!
! Global Reanalysis on T62 Gaussian grid is 192x94x28
!
!-----------------------------------------------------------------------

  INTEGER :: glreannx,glreanny,glreannz
  PARAMETER ( glreannx = 192, glreanny = 94, glreannz = 28 )

!-----------------------------------------------------------------------
!
! ARPS default for may20, 67x67x35
! NOTE:  ARPS dimensions are passed to CVT2VERIF by the input file. 
!
!-----------------------------------------------------------------------
!
!  parameter ( nx_ext = 67, ny_ext = 67, nz_ext = 35 )
!
!-----------------------------------------------------------------------
!
! Native RUC2 GRIB file (Grid #236) from OSO is 151x113x40
!
!-----------------------------------------------------------------------

  INTEGER :: rucn236nx,rucn236ny,rucn236nz
  PARAMETER (rucn236nx=151, rucn236ny=113, rucn236nz=40)

!-----------------------------------------------------------------------
!
! Isobaric RUC2 GRIB file (Grid #236) from OSO is 151x113x37
!
!-----------------------------------------------------------------------

  INTEGER :: rucp236nx,rucp236ny,rucp236nz
  PARAMETER (rucp236nx=151, rucp236ny=113, rucp236nz=37)

!-----------------------------------------------------------------------
!
! End MODULE extdims
!
!-----------------------------------------------------------------------

END MODULE extdims2
