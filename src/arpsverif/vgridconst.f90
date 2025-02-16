!
!##################################################################
!##################################################################
!######                                                      ######
!######                    VGRIDCONST.F90                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Include file for cvt2verif.f90, containing map projection information 
! on different verification grids.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp
! 11/2/1999
!
!-----------------------------------------------------------------------
!
MODULE vgridconst

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
!
! Map projection parameters for Eta GRIB #212
!
!-----------------------------------------------------------------------
!
  REAL,PARAMETER :: eta212dx = 40635.
  REAL,PARAMETER :: eta212dy = 40635.

  REAL,PARAMETER :: eta212ctrlat = 40.24
  REAL,PARAMETER :: eta212ctrlon = -101.00

  REAL,PARAMETER :: eta212latsw = 12.190
  REAL,PARAMETER :: eta212lonsw = -133.459

  INTEGER,PARAMETER :: eta212iproj = 2
  
  REAL,PARAMETER :: eta212scale = 1.

  REAL,PARAMETER :: eta212trlon = -95.0

  REAL,PARAMETER :: eta212latnot1 = 25.0
  REAL,PARAMETER :: eta212latnot2 = 25.0

!-----------------------------------------------------------------------
!
! Map projection parameters for RUC2 GRIB #236
!
!-----------------------------------------------------------------------
!

  REAL,PARAMETER :: ruc236dx = 40635.
  REAL,PARAMETER :: ruc236dy = 40635.

  REAL,PARAMETER :: ruc236ctrlat = 39.60
  REAL,PARAMETER :: ruc236ctrlon = -98.67

  REAL,PARAMETER :: ruc236latsw = 16.2810
  REAL,PARAMETER :: ruc236lonsw = -126.1387

  INTEGER,PARAMETER :: ruc236iproj = 2

  REAL,PARAMETER :: ruc236scale = 1.

  REAL,PARAMETER :: ruc236trlon = -95.

  REAL,PARAMETER :: ruc236latnot1 = 25.0
  REAL,PARAMETER :: ruc236latnot2 = 25.0

END MODULE vgridconst
