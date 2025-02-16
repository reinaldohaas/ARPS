!
!  Lookup table parameters
!  Variables beginning with "f" are real versions of the
!  integer variables, they must match the value of the
!  corresponding integer value.
!  This is for the benefit of the Cray, which doesn't
!  allow a reference to function "float" in a parameter statement.
!
  INTEGER :: nrnglut,nazmlut,nelvlut
  REAL :: fnrnglut,fnazmlut
  PARAMETER (nrnglut=250,                                               &
            fnrnglut=250.,                                              &
             nazmlut=360,                                               &
            fnazmlut=360.,                                              &
             nelvlut=6)
!             nelvlut=15)
  REAL :: farange,topelv,rnglutint,azmlutint
  PARAMETER (farange=250000.,                                           &
             topelv=20.,                                                &
             rnglutint=(farange/fnrnglut),                              &
             azmlutint=(360./fnazmlut))
!
!  List of elevations
!  This is list is appropriate for the WSR-88D programmed
!  scan stretegies.
!
  REAL :: elvlut(nelvlut)
!  DATA elvlut /0.0, 0.5, 1.5, 2.5, 3.5, 4.5,5.5,6.5,7.5,                &
!               8.5,10.0,12.0,14.0,16.0,19.5/
  DATA elvlut /0.0, 0.5, 1.5, 2.5, 3.5, 4.5/
!
  REAL :: radri(0:nrnglut,0:nazmlut,nelvlut)
  REAL :: radrj(0:nrnglut,0:nazmlut,nelvlut)
  REAL :: radrk(0:nrnglut,0:nazmlut,nelvlut)
  COMMON /lookup_rad/ radri,radrj,radrk
!
!  Lookup arrays at gridpoints
!
  REAL :: gridazm(nx,ny)
  REAL :: gridrng(nx,ny)
