!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                      RADCST.INC                      ######
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
!  Include file 'radcst.inc' for ARPS
!
!  This file initializes values of some constants used for radiation
!  calculation that remain unchanged throughout a model run
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/18/1996
!
!  MODIFICATION HISTORY:
!
!  08/14/1996 (Yuhe Liu)
!    Combined radims.inc and radcst.inc into radcst.inc. Defined the
!    radiation working arrays in the ARPS main program
!
!  04/09/1997 (Yuhe Liu)
!    Changed the indeces of working arrays from parameters to
!    variables and listed them in common blocks. This is to cooperate
!    with nested runs which have different dimension sizes for
!    different grids.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Define buffer array for option radopt=2. For option radopt=1, the
!  buffer size should be set to 1 to reduce the memory allocation.
!  Otherwise, the dimension sizes should be the same as nx, ny, and
!  nz, and the buffer size should be larger than the total size of 27
!  2-d arrays and 44 3-d arrays.
!
!    n2d_radiat  number of 2-d arrays in the buffer
!    n3d_radiat  number of 3-d arrays in the buffer
!
!  The 2-d arrays should be always at the beginning of radbuf and
!  the 3-d then follow.
!
!-----------------------------------------------------------------------
!
  INTEGER :: n2d_radiat  ! number of 2-d arrays in the buffer
  INTEGER :: n3d_radiat  ! number of 3-d arrays in the buffer

  PARAMETER ( n2d_radiat = 27 )
  PARAMETER ( n3d_radiat = 44 )
!
!-----------------------------------------------------------------------
!
!  Define the indeces that specify the positions of 2-D arrays in a
!  buffer array which carrys some of the radiation variables
!  calculated and/or used for both radiation options, radopt=1, or 2.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nrad2d   ! Number of 2-D arrays in radiation buffer array
  PARAMETER ( nrad2d = 10 )
           ! Buffur array to carry the variables calculated or used in
! radiation calculation. The last index defines the variables:
! 1  = nrsirbm,  Solar IR surface albedo for beam
! 2  = nrsirdf,  Solar IR surface albedo for diffuse
! 3  = nrsuvbm,  Solar UV surface albedo for beam
! 4  = nrsuvdf,  Solar UV surface albedo for diffuse
! 5  = ncosz,    Cosine of zenith
! 6  = ncosss,   Cosine of angle between sun light and
!                terrain slope
! 7  = nfdirir,  all-sky direct downward IR flux
!                (0.7-10 micron) at the surface
! 8  = nfdifir,  all-sky diffuse downward IR flux
           !                at the surface
! 9  = nfdirpar, all-sky direct downward par flux
!                (0.4-0.7 micron) at the surface
! 10 = nfdifpar, all-sky diffuse downward par flux
!                at the surface

  INTEGER :: nrsirbm  ! Index for solar IR surface albedo for beam
  PARAMETER ( nrsirbm  = 1 )

  INTEGER :: nrsirdf  ! Index for solar IR surface albedo for diffuse
  PARAMETER ( nrsirdf  = 2 )

  INTEGER :: nrsuvbm  ! Index for solar UV surface albedo for beam
  PARAMETER ( nrsuvbm  = 3 )

  INTEGER :: nrsuvdf  ! Index for solar UV surface albedo for diffuse
  PARAMETER ( nrsuvdf  = 4 )

  INTEGER :: ncosz    ! Index for cosine of solar zenith angle
  PARAMETER ( ncosz = 5 )

  INTEGER :: ncosss   ! Index for cosine of angle between sun light and
                      ! terrain slope
  PARAMETER ( ncosss = 6 )

  INTEGER :: nfdirir  ! Index for all-sky direct downward IR flux
                      !           (0.7-10 micron) at the surface
  PARAMETER ( nfdirir = 7 )

  INTEGER :: nfdifir  ! Index for all-sky diffuse downward IR flux
                      !           at the surface
  PARAMETER ( nfdifir = 8 )

  INTEGER :: nfdirpar ! Index for all-sky direct downward and par flux
                      !           (0.4-0.7 micron) at the surface
  PARAMETER ( nfdirpar = 9 )

  INTEGER :: nfdifpar ! Index for all-sky diffuse downward and par flux
                      !           at the surface
  PARAMETER ( nfdifpar = 10 )
!
!-----------------------------------------------------------------------
!
!  Define co2 concentration in mppv
!
!-----------------------------------------------------------------------
!
  REAL :: co2

  PARAMETER ( co2 = 300.0E-6 )  ! in mppv
!
!-----------------------------------------------------------------------
!
!  Define parameters for solid water components
!
!-----------------------------------------------------------------------
!
  REAL :: roqr,tnw,roqs,tns,roqg,tng

  PARAMETER ( roqr = 1.00   )
  PARAMETER ( tnw  = 0.08   )
  PARAMETER ( roqs = 0.10   )
  PARAMETER ( tns  = 0.03   )
  PARAMETER ( roqg = 0.913  )
  PARAMETER ( tng  = 0.0004 )
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine IRRAD and the subroutine IRRAD calls.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ir2d1,ir2d2,ir2d3,ir2d4,ir2d5,ir2d6     ! for 2-d arrays

  COMMON /ir2dcmn/ ir2d1,ir2d2,ir2d3,ir2d4,ir2d5,ir2d6

  INTEGER :: ir3d1, ir3d2, ir3d3, ir3d4, ir3d5     ! for 3-d arrays
  INTEGER :: ir3d6, ir3d7, ir3d8, ir3d9, ir3d10    ! for 3-d arrays
  INTEGER :: ir3d11,ir3d12,ir3d13,ir3d14,ir3d15    ! for 3-d arrays
  INTEGER :: ir3d16

  COMMON /ir3dcmn/ ir3d1, ir3d2, ir3d3, ir3d4, ir3d5,                   &
                   ir3d6, ir3d7, ir3d8, ir3d9, ir3d10,                  &
                   ir3d11,ir3d12,ir3d13,ir3d14,ir3d15,                  &
                   ir3d16

  INTEGER :: ir4d1, ir4d2                          ! for 4-d arrays

  COMMON /ir4dcmn/ ir4d1, ir4d2

  INTEGER :: ir5d1                                 ! for 5-d arrays

  COMMON /ir5dcmn/ ir5d1
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine SOLIR, SOLUV, and CLDFLX.
!
!-----------------------------------------------------------------------
!
  INTEGER :: so2d1, so2d2, so2d3, so2d4, so2d5     ! for 2-d arrays
  INTEGER :: so2d6, so2d7, so2d8, so2d9, so2d10
  INTEGER :: so2d11,so2d12,so2d13,so2d14,so2d15
  INTEGER :: so2d16,so2d17,so2d18,so2d19,so2d20
  INTEGER :: so2d21,so2d22,so2d23

  COMMON /so2dcmn/ so2d1, so2d2, so2d3, so2d4, so2d5,                   &
                   so2d6, so2d7, so2d8, so2d9, so2d10,                  &
                   so2d11,so2d12,so2d13,so2d14,so2d15,                  &
                   so2d16,so2d17,so2d18,so2d19,so2d20,                  &
                   so2d21,so2d22,so2d23

  INTEGER :: so3d1, so3d2, so3d3, so3d4, so3d5     ! for 3-d arrays
  INTEGER :: so3d6, so3d7, so3d8, so3d9, so3d10    ! for 3-d arrays
  INTEGER :: so3d11,so3d12,so3d13,so3d14           ! for 3-d arrays

  COMMON /so3dcmn/ so3d1, so3d2, so3d3, so3d4, so3d5,                   &
                   so3d6, so3d7, so3d8, so3d9, so3d10,                  &
                   so3d11,so3d12,so3d13,so3d14

  INTEGER :: so4d1, so4d2, so4d3, so4d4, so4d5     ! for 4-d arrays

  COMMON /so4dcmn/ so4d1, so4d2, so4d3, so4d4, so4d5

  INTEGER :: so5d1, so5d2, so5d3, so5d4, so5d5     ! for 5-d arrays

  COMMON /so5dcmn/ so5d1, so5d2, so5d3, so5d4, so5d5
