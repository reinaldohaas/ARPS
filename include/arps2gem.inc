!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                     ARPS2GEM.INC                     ######
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
!  PURPOSE:
!
!  Include file 'arps2gem.inc' for the arps2gem conversion program
!  from arps data to GEMPAK formatted data.
!
!  This file defines the number of pressure levels and
!  the actual pressure levels for which to interpolate in arps2gem
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jonathan Case
!  3/19/98
!
!-----------------------------------------------------------------------
!
!
!  nprgem:  number of pressure levels for GEMPAK file.
!
!  iprgem:  actual pressure levels defined by the user (mb)
!
  INTEGER :: nprgem
!
  PARAMETER (nprgem = 19)
!
  INTEGER :: iprgem(nprgem)
!
  DATA iprgem /1000, 950, 900, 850, 800, 750, 700,                      &
                650, 600, 550, 500, 450, 400, 350,                      &
                300, 250, 200, 150, 100/

!
!  end of include file
!
