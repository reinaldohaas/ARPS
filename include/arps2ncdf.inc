!
!##################################################################
!##################################################################
!######                                                      ######
!######                     ARPS2CDF.INC                     ######
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
!  Include file 'arps2cdf.inc' for the arps2gem conversion program
!  from arps data to GEMPAK formatted data.
!
!  This file defines the number of pressure levels and
!  the actual pressure levels for which to interpolate in arps2gem
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, based on arps2gem by John Case
!  2/15/99
!
!-----------------------------------------------------------------------
!
!
!  nprgem:  number of pressure levels for netCDF file.
!
!  iprlvl:  pressure levels defined by the user (mb)
!
  INTEGER :: nprlvl
!
  PARAMETER (nprlvl = 19)
!
  INTEGER :: iprlvl(nprlvl)
!
  DATA iprlvl /1000, 950, 900, 850, 800, 750, 700,                      &
                650, 600, 550, 500, 450, 400, 350,                      &
                300, 250, 200, 150, 100/
!
!-----------------------------------------------------------------------
!
!  End of include file
!
!-----------------------------------------------------------------------
!
