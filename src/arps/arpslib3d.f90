!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RANARY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ranary(nx,ny,nx1,nx2,ny1,ny2,iseed,amplit, rantp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate a 2-D array of machine-independent random numbers
!  between -amplit and +amplit with average value equal to zero.
!  The input parameter, iseed, must be a negative integer number
!  which needs to be defined only once in the calling subroutine
!  program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V.Wong and X.Song
!  7/20/1992.
!
!  Modification history:
!  7/25/1992. (MX)
!  Modified to allow assignment on a portion of the array.
!
!  12/11/1992 (MX)
!  Bug fix to the index bounds in loop 10.
!
!  9/1/94 (J. Levit & Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx         Number of random numbers in the x-direction
!    ny         Number of random numbers in the y-direction
!    nx1,nx2    Begin and end 1st array indicies specifying a
!               subdomain in which the average is set to zero
!    ny1,ny2    Begin and end 2nd array indicies specifying a
!               subdomain in which the average is set to zero
!    iseed      an arbitrary negative integer as a seed for a
!               sequence of random numbers
!    amplit     The generated numbers stay within the range
!               [-amplit,amplit]
!
!  OUTPUT:
!
!    rantp      The array for storing random numbers
!
!  WORK ARRAY:
!
!    tem1      Temporary working array
!
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations

  INCLUDE 'mp.inc'

  INTEGER :: i,j
  INTEGER :: nx,ny        ! Number of random numbers in each direction
  INTEGER :: nx1,nx2      ! Begin and end 1st array indecies specifying
                          ! a subdomain in which the average is set to
                          ! zero
  INTEGER :: ny1,ny2      ! Begin and end 2nd array indecies specifying
                          ! a subdomain in which the average is set to
                          ! zero
  INTEGER :: iseed        ! The seed for random number generation
  REAL    :: amplit       ! The generated numbers stay within
                          ! [-amplit, amplit]
  REAL    :: rantp (nx,ny)   ! Output array for storing the numbers
  REAL    :: ran3            ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: wv,ranave,rantpmax

  REAL, ALLOCATABLE :: tem1 (:,:)   ! Temporary working array for
                                    ! storing the numbers
  INTEGER           :: nxlg, nylg
  INTEGER           :: istatus
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

   nxlg = (nx-3)*nproc_x+3    ! For serial runs, this is just (nx,ny)
   nylg = (ny-3)*nproc_y+3

   ALLOCATE(tem1(nxlg,nylg), STAT = istatus)
   CALL check_alloc_status(istatus, "ranary:tem1")

!  Initialize tem1 to zero
   tem1 = 0.0
!
!-----------------------------------------------------------------------
!
!  Generate random numbers between 0 and 1, for a given iseed.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    ranave=0.0
    DO j=ny1,ny2
      DO i=nx1,nx2
        tem1(i,j)=ran3(iseed)
        ranave=ranave+tem1(i,j)
      END DO
    END DO
    ranave=ranave/(FLOAT( (nx2-nx1+1)*(ny2-ny1+1) ))
!
!-----------------------------------------------------------------------
!
!  Adjust the random numbers so that the new plane average value
!  equals zero and the random numbers stay within the range
!  [-amplit, amplit].
!
!-----------------------------------------------------------------------
!
    rantpmax=0.0
    DO j=ny1,ny2
      DO i=nx1,nx2
        rantpmax=AMAX1(rantpmax,ABS(ranave-tem1(i,j)))
      END DO
    END DO

    wv=amplit/rantpmax

    DO j=ny1,ny2
      DO i=nx1,nx2
        tem1(i,j)=(tem1(i,j)-ranave)*wv
      END DO
    END DO

  ENDIF !end of processor 0 section

  CALL mpisplit2d(tem1,nx,ny,rantp)

  DEALLOCATE(tem1)

  RETURN
END SUBROUTINE ranary

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION RAN3                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

REAL FUNCTION ran3(iseed)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generates a random number between 0 and 1 by feeding
!  a negative integer iseed.
!
!  Reference: "Seminumerical Algorithms" by Donald Knuth
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    iseed      an arbitrary negative integer as a seed for a
!               sequence of random numbers
!
!  OUTPUT:
!
!    ran3       A random number between 0 and 1.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: iseed        ! The seed for random number generation
!  REAL    :: ran3        ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: mbig,mseed,mz,k,ii,inext,inextp,i,iff,mk,mj
  REAL :: fac

  PARAMETER (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  INTEGER :: ma(55)
  SAVE iff, inext, inextp, ma
  DATA iff /0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!----------------------------------------------------------------------
!
!  Initialize the sequence of random numbers between 0 and 1,
!  using iseed.
!
!----------------------------------------------------------------------
!
  IF (iseed < 0.OR.iff == 0) THEN
    iff=1
    mj=mseed-IABS(iseed)
    mj=MOD(mj,mbig)
    ma(55)=mj
    mk=1
    DO i=1,54
      ii=MOD(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      IF (mk < mz) mk=mk+mbig
      mj=ma(ii)
    END DO

    DO k=1,4
      DO i=1,55
        ma(i)=ma(i)-ma(1+MOD(i+30,55))
        IF (ma(i) < mz) ma(i)=ma(i)+mbig
      END DO
    END DO

    inext=0
    inextp=31
    iseed=1
  END IF
!
!----------------------------------------------------------------------
!
!  Start to generate a random number.
!
!----------------------------------------------------------------------
!
  inext=inext+1
  IF (inext == 56) inext=1
  inextp=inextp+1
  IF (inextp == 56) inextp=1
  mj=ma(inext)-ma(inextp)
  IF (mj < mz) mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac

  RETURN
  END FUNCTION ran3
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SWAP4BYTE                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE swap4byte(a,n)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reverse order of bytes in integer*4 words
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: n
  INTEGER*4   a(n)

  INTEGER :: i

  INTEGER*4   itema
  CHARACTER (LEN=1) :: ctema(4)
  CHARACTER (LEN=1) :: ctemb

  EQUIVALENCE ( itema,ctema(1) )
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i = 1,n
    itema    = a(i)
    ctemb    = ctema(4)
    ctema(4) = ctema(1)
    ctema(1) = ctemb
    ctemb    = ctema(3)
    ctema(3) = ctema(2)
    ctema(2) = ctemb
    a(i)     = itema
  END DO

  RETURN
END SUBROUTINE swap4byte
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE FINDVARTYPE                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE findvartype(iendn,itypec,lw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find the type of endian, type of character set and length of
!  machine word.  Renamed version of subroutine w3fi04 from nmcdecode.f90.
!
!-----------------------------------------------------------------------
!
!SUBROUTINE w3fi04(iendn,itypec,lw)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3FI04         FIND WORD SIZE, ENDIAN, CHARACTER SET
!   PRGMNR: JONES,R.E.       ORG: W/NMC42       DATE: 94-10-07
!
! ABSTRACT: SUBROUTINE COMPUTES WORD SIZE, THE TYPE OF CHARACTER
!   SET, ASCII OR EBCDIC, AND IF THE COMPUTER IS BIG-ENDIAN, OR
!   LITTLE-ENDIAN.
!
! PROGRAM HISTORY LOG:
!   94-10-07  R.E.JONES
!
! USAGE:  CALL W3FI04 (IENDN, ITYPEC, LW)
!
!   OUTPUT ARGUMENT LIST:
!  IENDN     -  INTEGER FOR BIG-ENDIAN OR LITTLE-ENDIAN
!               = 0   BIG-ENDIAN
!               = 1   LITTLE-ENDIAN
!               = 2   CANNOT COMPUTE
!  ITYPEC    -  INTEGER FOR TYPE OF CHARACTER SET
!               = 0   ASCII  CHARACTER SET
!               = 1   EBCDIC CHARACTER SET
!               = 2   NOT ASCII OR EBCDIC
!  LW        -  INTEGER FOR WORDS SIZE OF COMPUTER IN BYTES
!               = 4   FOR 32 BIT COMPUTERS
!               = 8   FOR 64 BIT COMPUTERS
!
! ATTRIBUTES:
!   LANGUAGE: SiliconGraphics 3.5 FORTRAN 77
!   MACHINE:  SiliconGraphics IRIS-4D/25, 35, INDIGO, INDY
!
!$$$
!
  INTEGER :: itest1
  INTEGER :: itest2
  INTEGER :: itest3
  INTEGER :: iendn
  INTEGER :: itypec
  INTEGER :: lw
!
  CHARACTER (LEN=8) :: ctest1
  CHARACTER (LEN=8) :: ctest2
  CHARACTER (LEN=1) :: ctest3(8)
  CHARACTER (LEN=1) :: BLANK
!
  EQUIVALENCE   (ctest1,itest1),(ctest2,itest2)
!
  EQUIVALENCE   (itest3,ctest3(1))
!
  DATA  ctest1/'12345678'/
  DATA  itest3/z'01020304'/
  DATA  BLANK /' '/
!
  SAVE
!
!  TEST FOR TYPE OF CHARACTER SET
!  BLANK IS 32 (20 HEX) IN ASCII, 64 (40 HEX) IN EBCDEC
!
  IF (ICHAR(BLANK) == 32) THEN
    itypec = 0
  ELSE IF (ICHAR(BLANK) == 64) THEN
!
!  COMPUTER IS PROBABLY AN IBM360, 370, OR 390 WITH
!  A 32 BIT WORD SIZE, AND BIG-ENDIAN.
!
    itypec = 1
  ELSE
    itypec = 2
  END IF
!
!  TEST FOR WORD SIZE, SET LW TO 4 FOR 32 BIT COMPUTER,
!  8 FOR FOR 64 BIT COMPUTERS
!
  itest2 = itest1
  IF (ctest1 == ctest2) THEN
!
!  COMPUTER MAY BE A CRAY, OR COULD BE DEC VAX ALPHA
!  OR SGI WITH R4000, R4400, R8800 AFTER THEY CHANGE
!  FORTRAN COMPILERS FOR 64 BIT INTEGER.
!
    lw = 8
  ELSE
    lw = 4
  END IF
!
!  IF CHARACTER SET IS NOT ASCII OR EBCDIC SET IENDN = 2
!  CAN NOT TEST FOR ENDIAN TYPE
!
  IF (itypec == 2) THEN
    iendn = 2
    RETURN
  END IF
!
!  USING ITEST3 WITH Z'01020304' EQUIVALNCED TO CTEST3
!  ON A 32 BIT BIG-ENDIAN COMPUTER 03 IS IN THE 3RD
!  BYTE OF A 4 BYTE WORD. ON A 32 BIT LITTLE-ENDIAN
!  COMPUTER IT IS IN 2ND BYTE.
!  ON A 64 BIT COMPUTER Z'01020304' IS RIGHT ADJUSTED IN
!  A 64 BIT WORD, 03 IS IN THE 7TH BYTE.  ON A LITTLE-
!  ENDIAN 64 BIT COMPUTER IT IS IN THE 2ND BYTE.
!
  IF (lw == 4) THEN
    IF (ICHAR(ctest3(3)) == 3) THEN
      iendn = 0
    ELSE IF (ICHAR(ctest3(3)) == 2) THEN
      iendn = 1
    ELSE
      iendn = 2
    END IF
  ELSE IF (lw == 8) THEN
    IF (ICHAR(ctest3(7)) == 3) THEN
      iendn = 0
    ELSE IF (ICHAR(ctest3(2)) == 3) THEN
      iendn = 1
    ELSE
      iendn = 2
    END IF
  ELSE
    iendn = 2
  END IF
!
  RETURN
END SUBROUTINE findvartype
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE FIX_LAKE_ELIV               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fix_lake_eliv(h,lat,lon,nx,ny)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill in Great Lakes.  The original data set corresponds to the
!  bottom of the lakes.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Richard Carpenter.
!  2000/02/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT/OUPUT :
!
!    h        Height of the grid point.
!
!  INPUT :
!
!    lat      Latitude of the grid point.
!    lon      Longitude of the grid point.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!
!  OUTPUT:
!
!  WORK ARRAYS:
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny            ! x- y-Dimensions.

  REAL :: h(nx,ny)             ! Height
  REAL :: lat(nx,ny)           ! Latitude
  REAL :: lon(nx,ny)           ! Longitude

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: nlake, maxlake
  PARAMETER (maxlake=12)
  REAL :: lakelat1(maxlake), lakelat2(maxlake),                         &
      lakelon1(maxlake), lakelon2(maxlake), lakeelev(maxlake)

  REAL :: latpnt,lonpnt
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  PRINT *, 'ARPSTRN: Great Lakes terrain'
!
!-----------------------------------------------------------------------
!
!  Define the lat/lon boxes which contain each lake
!
!-----------------------------------------------------------------------
!
  k = 0

! Lake Superior
  k = k + 1
  lakelat1(k) = 46.4
  lakelat2(k) = 49.0
  lakelon1(k) = -92.5
  lakelon2(k) = -84.1
  lakeelev(k) = 182.9

! Lake Michigan
  k = k + 1
  lakelat1(k) = 41.0
  lakelat2(k) = 46.3
  lakelon1(k) = -88.1
  lakelon2(k) = -84.7
  lakeelev(k) = 176.5

! Lake Huron
  k = k + 1
  lakelat1(k) = 43.0
  lakelat2(k) = 46.4
  lakelon1(k) = -84.7
  lakelon2(k) = -79.5
  lakeelev(k) = 176.5

! Lake St. Clair
  k = k + 1
  lakelat1(k) = 42.2
  lakelat2(k) = 42.8
  lakelon1(k) = -83.0
  lakelon2(k) = -82.2
  lakeelev(k) = 174.6

! Lake Erie
  k = k + 1
  lakelat1(k) = 41.1
  lakelat2(k) = 43.0
  lakelon1(k) = -83.7
  lakelon2(k) = -78.6
  lakeelev(k) = 173.7

! Lake Ontario
  k = k + 1
  lakelat1(k) = 43.0
  lakelat2(k) = 45.6
  lakelon1(k) = -80.0
  lakelon2(k) = -76.0
  lakeelev(k) = 74.1

  nlake = k

  PRINT *, 'ARPSTRN: Number of Great Lakes: ', nlake
!
!-----------------------------------------------------------------------
!
!  Do no allow the elevation in a box to go below the lake's elevation.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      latpnt = lat(i,j)
      lonpnt = lon(i,j)
      IF(lonpnt > 180) lonpnt = lonpnt - 360

      DO k=1,nlake
        IF (latpnt > lakelat1(k) .AND. latpnt < lakelat2(k) .AND.       &
                lonpnt > lakelon1(k) .AND. lonpnt < lakelon2(k) .AND.   &
                h(i,j) < lakeelev(k) ) THEN
!        write (*,'(a,2i4,2f8.1,i4,2f7.1)')
!    :           "FIXELEVL",i,j,lakeelev(k),h(i,j),k,latpnt,lonpnt
          h(i,j) = lakeelev(k)
        END IF
      END DO
!    write (*,'(a,2i4,f8.1,12x,2f7.1)')
!    :           "FIXELEV ",i,j,h(i,j),latpnt,lonpnt

    END DO
  END DO


  RETURN
END SUBROUTINE fix_lake_eliv
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE EXTMNSND                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extmnsnd(nx,ny,nx_ext,ny_ext,nz_ext,lvlprof,ovrlay,          &
                    iextmn,iextmx,jextmn,jextmx,kextmn,kextmx,          &
                    xscl,yscl,x_ext,y_ext,                              &
                    hgt_ext,pln_ext,t_ext,                              &
                    rhs_ext,u_ext,v_ext,                                &
                    zsnd,plsnd,psnd,tsnd,ptsnd,rhssnd,qvsnd,            &
                    rhosnd,usnd,vsnd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Finds a mean sounding of all external data within the ARPS
!  domain.  This is used to define the model base state.
!  The data are interpolated to constant height levels and
!  averaged in the horizontal. The final sounding is hydrostatic.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  August, 1995  Replaces MNSOUND
!
!  MODIFICATION HISTORY:
!  12/20/95  Keith Brewster Added code to insure at least one
!            point goes into creating the mean sounding.
!
!  2000/04/05 (Gene Bassett)
!  Added qvsnd array.
!
!  04/02/2001 (K. Brewster)
!  Corrected calculations for qvsnd array that were also impacting
!  rhssnd output array.
!
!  Also, added parameter rhmin=0.05 (5 percent) as a minimum relative
!  humidity in place of hardcoded check value of 1.0E-04.
!
!  20 June 2002 (Eric Kemp)
!  Bug fix for kext loop.
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
!  Sizing variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny
  INTEGER :: nx_ext,ny_ext,nz_ext
  INTEGER :: lvlprof

  INTEGER, INTENT(IN) :: ovrlay     ! Domain overlay type, for wrf2arps_mpi & ext2arps_mpi
                                    ! = 1, External domain may be overlay, ARPS domain
                                    !      should not have overlay, such as ext2arps_mpi
                                    ! = 2, ARPS domain can be overlay, but external domain
                                    !      must not overlay, such as wrf2arps_mpi
!
!-----------------------------------------------------------------------
!
!  Limits of calculations
!
!-----------------------------------------------------------------------
!
  INTEGER :: iextmn,iextmx
  INTEGER :: jextmn,jextmx
  INTEGER :: kextmn,kextmx
!
!-----------------------------------------------------------------------
!
!  ARPS grid variables
!
!-----------------------------------------------------------------------
!
  REAL :: xscl(nx)
  REAL :: yscl(ny)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext,ny_ext)
  REAL :: y_ext(nx_ext,ny_ext)

  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)         ! Height MSL
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)           ! Temperature (K)
  REAL :: pln_ext(nx_ext,ny_ext,nz_ext)         ! ln of total pressure
  REAL :: rhs_ext(nx_ext,ny_ext,nz_ext)         ! RHstar=SQRT(1.-RH)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)           ! u wind component
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)           ! v wind component
!
!-----------------------------------------------------------------------
!
!  Input sounding variables
!
!-----------------------------------------------------------------------
!
  REAL :: zsnd(lvlprof)
!
!-----------------------------------------------------------------------
!
!  Output sounding variables
!
!-----------------------------------------------------------------------
!
  REAL :: plsnd(lvlprof)         ! ln of pressure
  REAL :: psnd(lvlprof)          ! Pressure (Pa)
  REAL :: tsnd(lvlprof)          ! Temperature (K)
  REAL :: ptsnd(lvlprof)         ! Potential temperature (K)
  REAL :: rhssnd(lvlprof)        ! RHstar
  REAL :: qvsnd(lvlprof)         ! Specific humidity (kg/kg)
  REAL :: rhosnd(lvlprof)        ! Density
  REAL :: usnd(lvlprof)          ! u wind component
  REAL :: vsnd(lvlprof)          ! v wind component

!
!-----------------------------------------------------------------------
!
!  Output sounding variables, double precision version
!
!  Running a sum on lots of small numbers ends up with a precision
!  problem.  MPI and non-MPI runs will NEVER have the same answers.
!  Non-MPI runs will not even have the same answers in the numbers are
!  computed in a different order!
!
!-----------------------------------------------------------------------
!
  DOUBLE PRECISION :: ptsnd_dbl(lvlprof)
  DOUBLE PRECISION :: plsnd_dbl(lvlprof)
  DOUBLE PRECISION :: tsnd_dbl(lvlprof)
  DOUBLE PRECISION :: rhssnd_dbl(lvlprof)
  DOUBLE PRECISION :: usnd_dbl(lvlprof)
  DOUBLE PRECISION :: vsnd_dbl(lvlprof)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ibgn, jbgn, iend, jend
  REAL    :: xscl_min_lg, xscl_max_lg
  REAL    :: yscl_min_lg, yscl_max_lg
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Lapse rate and a constant for hydrostatic integration
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: gamma  = 0.0068         ! 6.8 degrees per km
  REAL, PARAMETER :: pconst = (g/rd)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: rhmin=0.05
  REAL, PARAMETER :: rhmax=1.0

  INTEGER :: i,j,k,kext,imid,jmid
  INTEGER :: knt,knttot,ktop,kbot
  REAL    :: wlow,whigh,tvbar,xmid,ymid,dist2,distmn,accept
  REAL    :: c1,c2,pres,qvsat,qvbot,qvtop,rh,epsilon

!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  epsilon=0.01*(zsnd(2)-zsnd(1))
!
!-----------------------------------------------------------------------
!
!  Zero-out all sounding arrays
!  ptsnd array is used temporarily to store a count.
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    ptsnd(k) = 0.
    plsnd(k) = 0.
    tsnd(k)  = 0.
    rhssnd(k)= 0.
    qvsnd(k) = 0.
    usnd(k)  = 0.
    vsnd(k)  = 0.

    ptsnd_dbl(k)  = 0.
    plsnd_dbl(k)  = 0.
    tsnd_dbl(k)   = 0.
    rhssnd_dbl(k) = 0.
    usnd_dbl(k)   = 0.
    vsnd_dbl(k)   = 0.
  END DO

!
!-----------------------------------------------------------------------
!
!  Find the external data grid location that is closest
!  to the middle of the grid.
!  For the case where no external grid point lies within the
!  grid this grid column will be used as the mean sounding.
!
!-----------------------------------------------------------------------
!

  xscl_min_lg = xscl(1)
  xscl_max_lg = xscl(nx-1)
  yscl_min_lg = yscl(1)
  yscl_max_lg = yscl(ny-1)

  call mpmax0(xscl_max_lg, xscl_min_lg)
  call mpmax0(yscl_max_lg, yscl_min_lg)

  xmid=0.5*(xscl_min_lg + xscl_max_lg)
  ymid=0.5*(yscl_min_lg + yscl_max_lg)

  distmn=(x_ext(1,1)-xmid)*(x_ext(1,1)-xmid) +                          &
         (y_ext(1,1)-ymid)*(y_ext(1,1)-ymid)
  imid=1
  jmid=1
  DO j=1,ny_ext
    DO i=1,nx_ext
      dist2=(x_ext(i,j)-xmid)*(x_ext(i,j)-xmid) +                       &
            (y_ext(i,j)-ymid)*(y_ext(i,j)-ymid)
      IF(dist2 < distmn) THEN
        distmn=dist2
        imid=i
        jmid=j
      END IF
    END DO
  END DO
  iextmn=MIN(iextmn,imid)
  iextmx=MAX(iextmx,imid)
  jextmn=MIN(jextmn,jmid)
  jextmx=MAX(jextmx,jmid)

!
!-----------------------------------------------------------------------
!
!  Look at all external points possibly within domain
!
!-----------------------------------------------------------------------
!
  knt=0
  knttot=0
!
!  Make sure we do the computations only once!
!
  ibgn = 1
  jbgn = 1
  iend = nx
  jend = ny
  IF (mp_opt > 0 .AND. ovrlay == 1) THEN  ! for ext2arps_mpi to work
    IF (loc_x > 1)         ibgn = 2
    IF (loc_x < nproc_x)   iend = nx-1
    IF (loc_y > 1)         jbgn = 2
    IF (loc_y < nproc_y)   jend = ny-1
!  ELSE                   ! for wrf2arps_mpi to work
!    ibgn = 1             ! Suppose iextmn, iextmx, jextmn, jextmx do not
!    jbgn = 1             ! have overlaps. Then ovelaps with ARPS dimensions
!    iend = nx            ! will do no harm.
!    jend = ny
  END IF

  DO j = jextmn,jextmx
    DO i = iextmn,iextmx
!
!-----------------------------------------------------------------------
!
!  Is this point within the ARPS domain?
!  Since the ARPS grid is Cartesian, need only compare
!  the external grid coordinates to the x and y limits
!  of the ARPS grid.
!
!-----------------------------------------------------------------------
!
      knttot=knttot+1
      IF( (x_ext(i,j) >= xscl(ibgn) .AND. x_ext(i,j) < xscl(iend)) .AND. &
          (y_ext(i,j) >= yscl(jbgn) .AND. y_ext(i,j) < yscl(jend)) )THEN
!
!  From above:
!	For the case where no external grid point lies with the
!	grid this grid column will be used as the mean sounding.
!
!  This really doesn't make sense to me.  If there is no overlap between
!  the grid and the external data, then nothing useful should be expected
!  to happen.  That part of the code is commented out.  If it didn't have an
!  MPI issue (matches in all processors), I'd leave it in.  KWT
!
!          y_ext(i,j) >= yscl(jsub) .AND. y_ext(i,j) <= yscl(ny) )  .OR.   &
!         (i == imid .AND. j == jmid )) THEN

        knt=knt+1
!
!-----------------------------------------------------------------------
!
!  Interpolate external data in vertical onto profile
!  arrays.
!
!-----------------------------------------------------------------------
!
        DO k=1,lvlprof
          IF(zsnd(k) >= hgt_ext(i,j,kextmx)) EXIT
          IF((hgt_ext(i,j,1)-zsnd(k)) < epsilon) THEN
            DO kext=MAX(kextmn,2),kextmx-1           ! EMK 20 June 2002
              IF(hgt_ext(i,j,kext) >= zsnd(k)) EXIT
            END DO
            whigh=(zsnd(k)-hgt_ext(i,j,kext-1))/                        &
                  (hgt_ext(i,j,kext)-hgt_ext(i,j,kext-1))
            wlow=1.-whigh
            ptsnd_dbl(k)=ptsnd_dbl(k)+1.
            plsnd_dbl(k)=plsnd_dbl(k)+                                  &
                  whigh*pln_ext(i,j,kext)+wlow*pln_ext(i,j,kext-1)
            tsnd_dbl(k)=tsnd_dbl(k)+                                    &
                  whigh*t_ext(i,j,kext) + wlow*t_ext(i,j,kext-1)
            rhssnd_dbl(k)=rhssnd_dbl(k)+                                &
                  whigh*rhs_ext(i,j,kext)+wlow*rhs_ext(i,j,kext-1)
            usnd_dbl(k)=usnd_dbl(k)+                                    &
                  whigh*u_ext(i,j,kext) + wlow*u_ext(i,j,kext-1)
            vsnd_dbl(k)=vsnd_dbl(k)+                                    &
                  whigh*v_ext(i,j,kext) + wlow*v_ext(i,j,kext-1)
          END IF
        END DO
      END IF
    END DO
  END DO

!
! mpi sum  ptsnd(k), knt, knttot
!          plsnd(k),tsnd(k),rhssnd(k),usnd(k),vsnd(k)
!
  IF(mp_opt > 0) THEN

    CALL mptotali(knt)
    CALL mptotali(knttot)

    CALL mpsumdp(ptsnd_dbl, lvlprof)
    CALL mpsumdp(plsnd_dbl, lvlprof)
    CALL mpsumdp(tsnd_dbl,  lvlprof)
    CALL mpsumdp(rhssnd_dbl,lvlprof)
    CALL mpsumdp(usnd_dbl,  lvlprof)
    CALL mpsumdp(vsnd_dbl,  lvlprof)
  END IF

!
!  Return to single precision.
!

    DO k=1,lvlprof
      ptsnd(k)  = ptsnd_dbl(k)
      plsnd(k)  = plsnd_dbl(k)
      tsnd(k)   = tsnd_dbl(k)
      rhssnd(k) = rhssnd_dbl(k)
      usnd(k)   = usnd_dbl(k)
      vsnd(k)   = vsnd_dbl(k)
    END DO

!
! NOTE: plsnd, tsnd, etc should be exactly the same as
!       non-mpi mode because of the problem variables were made
!       double precision, which get around the round-off error.
!       The variables summed had sufficient precision.
!
!       "knttot" will be different in mpi and non-mpi modes, however, this
!       is totally harmless, as the last use is the WRITE below.  The numbers
!       can be made to be identical, however, doing it "the right way" makes
!       the run time longer for something that isn't of any value.
!


  IF(myproc == 0) WRITE(6,'(/a,i10,a,/12x,a,i10,a/)')                   &
       '  extmnsnd found ',knt,' points within ARPS domain',            &
       '  of ',knttot,' points in external domain checked.'

!
!-----------------------------------------------------------------------
!
!  Find lowest height with data
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a)') '  Finding range of mean sounding data ...'
  accept=0.3*knt
  DO k=1,lvlprof-1
    IF(myproc == 0) WRITE(6,'(a,f10.2,2(a,f10.0))')                     &
              ' z = ',zsnd(k),' knt = ',ptsnd(k),' accept = ',accept
    IF(ptsnd(k) > accept) EXIT
  END DO
  kbot=k
!
!-----------------------------------------------------------------------
!
!  Find highest height with data
!
!-----------------------------------------------------------------------
!
  DO k=lvlprof,2,-1
    IF(myproc == 0) WRITE(6,'(a,f10.2,2(a,f10.0))')                     &
             ' z = ',zsnd(k),' knt = ',ptsnd(k),' accept = ',accept
    IF(ptsnd(k) > accept) EXIT
  END DO
  ktop=k

  IF(myproc == 0) WRITE(6,'(a,f10.2,a,f10.2,a)')                        &
               ' Height of external data for mean spans from ',         &
                 zsnd(kbot),' to ',zsnd(ktop),' meters.'
!
!-----------------------------------------------------------------------
!
!  Divide through to find average.  We return to single precision here.
!
!-----------------------------------------------------------------------
!
  DO k=kbot,ktop
    plsnd(k)  = plsnd(k)/ptsnd(k)
    tsnd(k)   = tsnd(k)/ptsnd(k)
    rhssnd(k) = rhssnd(k)/ptsnd(k)
    usnd(k)   = usnd(k)/ptsnd(k)
    vsnd(k)   = vsnd(k)/ptsnd(k)
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "below-ground"
!  Use a constant lapse rate, gamma.
!  plsnd is a sort of first-guess log(pressure), needed for qv
!  calculation from rhstar.
!
!-----------------------------------------------------------------------
!
  pres = EXP(plsnd(kbot))
  rh=MAX(rhmin,rhmax-(rhssnd(kbot)*rhssnd(kbot)))
  qvsat = f_qvsat( pres, tsnd(kbot) )
  qvbot=rh*qvsat
  c1=g/(rd*gamma)
  DO k=kbot-1,1,-1
    tsnd(k)=tsnd(kbot)-gamma*(zsnd(k)-zsnd(kbot))
    plsnd(k)=plsnd(kbot)+                                               &
        c1*ALOG((tsnd(kbot)-gamma*(zsnd(k)-zsnd(kbot)))/tsnd(kbot))
    psnd(k) = EXP(plsnd(k))
    qvsat = f_qvsat( psnd(k), tsnd(k) )
    rh=qvbot/qvsat
    rhssnd(k)=SQRT(MAX(0.,(rhmax-rh)))
    usnd(k)=usnd(kbot)
    vsnd(k)=vsnd(kbot)
  END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "above-top"
!  Use a constant temperature. We're assuming stratosphere here.
!
!-----------------------------------------------------------------------
!
  pres = EXP(plsnd(ktop))
  rh=MAX(rhmin,rhmax-(rhssnd(ktop)*rhssnd(ktop)))
  qvsat = f_qvsat( pres, tsnd(ktop) )
  qvtop=rh*qvsat
  c2=g/rd
  DO k=ktop+1,lvlprof
    tsnd(k)=tsnd(ktop)
    plsnd(k)=plsnd(ktop)-c2*(zsnd(k)-zsnd(ktop))/tsnd(ktop)
    psnd(k) = EXP(plsnd(k))
    qvsat = f_qvsat( psnd(k), tsnd(k) )
    rh=qvtop/qvsat
    rhssnd(k)=SQRT(MAX(0.,(rhmax-rh)))
    usnd(k)=usnd(ktop)
    vsnd(k)=vsnd(ktop)
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate qv profile from RH-star.
!  Temporarily use rhosnd to store virtual temperature
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    pres = EXP(plsnd(k))
    rh=MAX(rhmin,rhmax-(rhssnd(k)*rhssnd(k)))
    qvsat = f_qvsat( pres, tsnd(k) )
    qvsnd(k)=rh*qvsat
    rhosnd(k) = tsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Make sure that the sounding is hydrostatic by integrating
!  from kbot. rhosnd is really virtual temperature here.
!
!-----------------------------------------------------------------------
!
  psnd(kbot)=EXP(plsnd(kbot))
  DO k=kbot-1,1,-1
    tvbar=0.5*(rhosnd(k+1)+rhosnd(k))
    psnd(k)=psnd(k+1)*EXP(pconst*(zsnd(k+1)-zsnd(k))/tvbar)
  END DO

  DO k=kbot+1,lvlprof
    tvbar=0.5*(rhosnd(k-1)+rhosnd(k))
    psnd(k)=psnd(k-1)*EXP(pconst*(zsnd(k-1)-zsnd(k))/tvbar)
  END DO
!
!-----------------------------------------------------------------------
!
!  Derived variable calculations
!  Compute density from virtual temperature.
!  Compute potential temperature.
!  Compute log of pressure.
!
!-----------------------------------------------------------------------
!
  DO k=1,lvlprof
    rhosnd(k)=psnd(k)/(rd*rhosnd(k))
    ptsnd(k)=tsnd(k)*((p0/psnd(k))**rddcp)
!    pres=exp(plsnd(k))
!    print *,'psnd old, psnd new: ',pres,psnd(k),(psnd(k)-pres)
    plsnd(k)=ALOG(psnd(k))
  END DO
!
  RETURN
END SUBROUTINE extmnsnd

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ARPSSTOP                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE arpsstop(comment,outerr)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Close down any necessary options (e.g. MPI) and stop execution of
!  the program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/21
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    comment    Comment string.
!    outerr     Output error indicator (0-send to standard out,
!                  1-send to standard error)
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  CHARACTER (LEN=*) :: comment     ! comment
  INTEGER :: outerr       ! error code (0-standard out, 1-standard error)

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  WRITE(6,*)'EMK: arpsstop: outerr = ',outerr

  IF (outerr == 0) THEN
    WRITE (6,*) trim(comment)
  ELSE
    WRITE (0,*) trim(comment)
  END IF

  IF (mp_opt > 0) THEN
    CALL mpexit(outerr)
  ELSE
    CALL exit(outerr)
  END IF

  STOP
END SUBROUTINE arpsstop

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTCOMMENT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE wrtcomment(comment,outerr)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out a message to standard our or standard error.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/21
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    comment    Comment string.
!    outerr     Output error indicator (0-send to standard out,
!                  1-send to standard error)
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  CHARACTER (LEN=*) :: comment     ! comment
  INTEGER :: outerr       ! error code (0-standard out, 1-standard error)

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (outerr == 0) THEN
    IF (myproc == 0) WRITE (6,*) trim(comment)
  ELSE
    WRITE (0,*) trim(comment)
  END IF

  RETURN
END SUBROUTINE wrtcomment

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ACCT_INIT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE acct_init

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize cpu and wall clock timers.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  REAL :: f_cputime
  DOUBLE PRECISION :: f_walltime

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use_acct = 1

  CALL init_walltime

  wall_init = f_walltime()
  cpu_init  = dble(f_cputime())

  acct_cpu_time  = 0d0
  acct_wall_time = 0d0

  current_acct   = 0
  interrupt_acct = 0

  RETURN
END SUBROUTINE acct_init

!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SET_ACCT                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE set_acct(account)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Start cpu and wall timers for the specified account number after
!  adding the values of the current account to its totals.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    account    Integer pointer to the account to start timing.
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: account

  REAL :: f_cputime
  DOUBLE PRECISION :: f_walltime

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

  DOUBLE PRECISION :: cpu0, wall0

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (use_acct == 0) RETURN

  cpu0  = dble(f_cputime())
  wall0 =      f_walltime()

  IF (current_acct /= 0) THEN  ! close off current account

    acct_cpu_time(current_acct)  = acct_cpu_time(current_acct)          &
                                 + cpu0 - cpu_acct
    acct_wall_time(current_acct) = acct_wall_time(current_acct)         &
                                 + wall0 - wall_acct

    IF (interrupt_acct /= 0) THEN  ! reset the interrupt
      CALL acct_interrupt(interrupt_acct)
    END IF

  ENDIF

  current_acct = account

  cpu_acct = cpu0
  wall_acct = wall0

  RETURN
END SUBROUTINE set_acct

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ACCT_INTERRUPT             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE acct_interrupt(account)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interrupt accounting for the main timer to start tracting another
!  account.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    account    Integer pointer to the account to start timing.  A value
!               of zero will turn off all but to total cpu/wall timers.
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: account

  REAL :: f_cputime
  DOUBLE PRECISION :: f_walltime

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

  DOUBLE PRECISION :: cpu0, wall0

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (use_acct == 0) RETURN

  cpu0  = dble(f_cputime())
  wall0 =      f_walltime()

  IF (interrupt_acct /= 0) THEN  ! stop an existing interrupt
    CALL acct_stop_inter
  ENDIF

  interrupt_acct = account

  cpu_inter  = cpu0
  wall_inter = wall0

  RETURN
END SUBROUTINE acct_interrupt

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ACCT_STOP_INTER            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE acct_stop_inter

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Stop the interrupt timer and adjust the main timer accordingly.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  REAL :: f_cputime
  DOUBLE PRECISION :: f_walltime

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

  DOUBLE PRECISION :: cpu0, wall0

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (use_acct == 0) RETURN

  cpu0  = dble(f_cputime())
  wall0 =      f_walltime()

  IF (interrupt_acct /= 0) THEN
    acct_cpu_time(current_acct)  =   acct_cpu_time(current_acct)        &
                                   - (cpu0 - cpu_inter)
    acct_wall_time(current_acct) =   acct_wall_time(current_acct)       &
                                   - (wall0 - wall_inter)

    acct_cpu_time(interrupt_acct)  = acct_cpu_time(interrupt_acct)      &
                                   + (cpu0 - cpu_inter)
    acct_wall_time(interrupt_acct) = acct_wall_time(interrupt_acct)     &
                                   + (wall0 - wall_inter)
  END IF

  interrupt_acct = 0

  RETURN
END SUBROUTINE acct_stop_inter

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE ACCT_FINISH                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE acct_finish

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Stop cpu/wall timers and calculate totals.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

  REAL :: f_cputime
  DOUBLE PRECISION :: f_walltime

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

  DOUBLE PRECISION :: cpu0, wall0
  !REAL :: time
  !INTEGER :: acct

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (use_acct == 0) RETURN

  cpu0  = dble(f_cputime())
  wall0 =      f_walltime()

  total_cpu  = cpu0  - cpu_init
  total_wall = wall0 - wall_init

  IF (interrupt_acct /= 0) THEN  ! stop an existing interrupt
    CALL acct_stop_inter
  END IF

  IF (current_acct /= 0) THEN
    CALL set_acct(0)
  END IF

!-----------------------------------------------------------------------
!
!  Sum up the CPU times used by all processors if on a multi-processor
!  system.
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    CALL mpsumdp(total_wall, 1)
    !time = real(total_wall)
    !CALL mptotal(time)
    !total_wall = time

    CALL mpsumdp(total_cpu, 1)
    !time = real(total_cpu)
    !CALL mptotal(time)
    !total_cpu = time

    CALL mpsumdp(acct_cpu_time,  max_acct)
    CALL mpsumdp(acct_wall_time, max_acct)
    !DO acct = 1,max_acct
    !  time = real(acct_cpu_time(acct))
    !  CALL mptotal(time)
    !  acct_cpu_time(acct) = dble(time)
    !
    !  time = real(acct_wall_time(acct))
    !  CALL mptotal(time)
    !  acct_wall_time(acct) = dble(time)
    !END DO
  END IF

  RETURN
END SUBROUTINE acct_finish

!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE ACCT_REPORT_ARPS             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE acct_report_arps

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Output the timing statistics for the ARPS model.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/09/14
!
!  MODIFICATION HISTORY:
!
!  13 March 2002 (Eric Kemp)
!  Added time accounting for WRF BMJ scheme.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE        ! Force explicit declarations

!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------

  REAL :: tc, tw

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL :: wt

  CHARACTER(LEN=6), PARAMETER :: mscheme(0:12) = (/ 'None  ',           &
                       'warmra','linice','nemice','LFO   ',             &
                       'WSM6WR','WSM6GR','WSM6N0',                      &
                       'MYSM  ','MYDMf ','MYDMd ','MYTM  ','MY2MOM' /)

  CHARACTER(LEN=6), PARAMETER :: cumscheme(0:5) = (/                    &
                       'NO    ','KUO1  ','KUO2  ','KF    ',             &
                       'WRFBMJ','KF_ETA' /)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL flush(6)

  ! A barrier before and after printing stats avoids spurious output
  ! in the body of the stats output.
  CALL mpbarrier

  IF (myproc == 0) THEN

    IF( mp_opt == 1 ) then
      wt = 1.0/( nproc_x*nproc_y )
    ELSE
      wt = 1.0
    ENDIF

    IF(total_cpu  == 0d0) total_cpu =1d0
    IF(total_wall == 0d0) total_wall=1d0
    tc = 100.0/total_cpu
    tw = 100.0/total_wall

    !
    !  Print out CPU usage.
    !
    WRITE (6,'(/a/)') "ARPS CPU Summary:"

    IF (mp_opt > 0) THEN
      WRITE  (6,'(a,/,a)') &
      ' Process             CPU time                   WALL CLOCK time (mean)', &
      ' -----------------   ----------------------     ----------------------'
    ELSE
      WRITE  (6,'(a,/,a)') &
      ' Process             CPU time                   WALL CLOCK time',        &
      ' -----------------   ----------------------     ----------------------'
    END IF

    WRITE(6,'(25(1x,a,2x,1p,e13.6,a,2x,0p,f6.2,a,4x,1p,e13.6,a,2x,0p,f6.2,a,/))') &
    'Initialization  :',acct_cpu_time(init_acct),                      's', &
     acct_cpu_time(init_acct)*tc,  '%', acct_wall_time(init_acct)*wt,  's', &
     acct_wall_time(init_acct)*tw, '%',                                 &
    'Data output     :',acct_cpu_time(output_acct),                    's', &
     acct_cpu_time(output_acct)*tc,'%', acct_wall_time(output_acct)*wt,'s', &
     acct_wall_time(output_acct)*tw,'%',                                &
    'Wind advection  :',acct_cpu_time(advuvw_acct),                    's', &
     acct_cpu_time(advuvw_acct)*tc,'%', acct_wall_time(advuvw_acct)*wt,'s', &
     acct_wall_time(advuvw_acct)*tw,'%',                                &
    'Scalar advection:',acct_cpu_time(advs_acct),                      's', &
     acct_cpu_time(advs_acct)*tc,  '%', acct_wall_time(advs_acct)*wt,  's', &
     acct_wall_time(advs_acct)*tw, '%',                                 &
    'Coriolis force  :',acct_cpu_time(coriol_acct),                    's', &
     acct_cpu_time(coriol_acct)*tc,'%', acct_wall_time(coriol_acct)*wt,'s', &
     acct_wall_time(coriol_acct)*tw,'%',                                &
    'Buoyancy term   :',acct_cpu_time(buoy_acct),                      's', &
     acct_cpu_time(buoy_acct)*tc,  '%', acct_wall_time(buoy_acct)*wt,  's', &
     acct_wall_time(buoy_acct)*tw,'%',                                  &
    'Misc Large tstep:',acct_cpu_time(tinteg_acct),                    's', &
     acct_cpu_time(tinteg_acct)*tc,'%', acct_wall_time(tinteg_acct)*wt,'s', &
     acct_wall_time(tinteg_acct)*tw,'%',                                &
    'Small time steps:',acct_cpu_time(smlstp_acct),                    's', &
     acct_cpu_time(smlstp_acct)*tc,'%', acct_wall_time(smlstp_acct)*wt,'s', &
     acct_wall_time(smlstp_acct)*tw,'%',                                &
    'Radiation       :',acct_cpu_time(rad_acct),                       's', &
     acct_cpu_time(rad_acct)*tc,   '%', acct_wall_time(rad_acct)*wt,   's', &
     acct_wall_time(rad_acct)*tw, '%',                                  &
    'Soil model      :',acct_cpu_time(soil_acct),                      's', &
     acct_cpu_time(soil_acct)*tc,  '%', acct_wall_time(soil_acct)*wt,  's', &
     acct_wall_time(soil_acct)*tw,'%',                                  &
    'Surface physics :',acct_cpu_time(sfcphy_acct),                    's', &
     acct_cpu_time(sfcphy_acct)*tc,'%', acct_wall_time(sfcphy_acct)*wt,'s', &
     acct_wall_time(sfcphy_acct)*tw,'%',                                &
    'Turbulence      :',acct_cpu_time(tmix_acct),                      's', &
     acct_cpu_time(tmix_acct)*tc,  '%', acct_wall_time(tmix_acct)*wt,  's', &
     acct_wall_time(tmix_acct)*tw,'%',                                  &
    'Comput. mixing  :',acct_cpu_time(cmix_acct),                      's', &
     acct_cpu_time(cmix_acct)*tc,  '%', acct_wall_time(cmix_acct)*wt,  's', &
     acct_wall_time(cmix_acct)*tw,'%',                                  &
    'Rayleigh damping:',acct_cpu_time(raydmp_acct),                    's', &
     acct_cpu_time(raydmp_acct)*tc,'%', acct_wall_time(raydmp_acct)*wt,'s', &
     acct_wall_time(raydmp_acct)*tw,'%',                                &
    'TKE src terms   :',acct_cpu_time(tkesrc_acct),                    's', &
     acct_cpu_time(tkesrc_acct)*tc,'%', acct_wall_time(tkesrc_acct)*wt,'s', &
     acct_wall_time(tkesrc_acct)*tw,'%',                                &
    'Gridscale precp.:',acct_cpu_time(qpfgrd_acct),                    's', &
     acct_cpu_time(qpfgrd_acct)*tc,'%', acct_wall_time(qpfgrd_acct)*wt,'s', &
     acct_wall_time(qpfgrd_acct)*tw,'%',                                &
    'Cumulus ('//cumscheme(cnvctopt)//'):',acct_cpu_time(cum_acct),    's', &
     acct_cpu_time(cum_acct)*tc,'%', acct_wall_time(cum_acct)*wt,      's', &
     acct_wall_time(cum_acct)*tw,'%',                                   &
    'Microph ('//mscheme(mphyopt)//'):',acct_cpu_time(microph_acct),    's',&
     acct_cpu_time(microph_acct)*tc,'%',acct_wall_time(microph_acct)*wt,'s',&
     acct_wall_time(microph_acct)*tw,'%',                               &
    'Hydrometero fall:',acct_cpu_time(qfall_acct),                     's', &
     acct_cpu_time(qfall_acct)*tc, '%', acct_wall_time(qfall_acct)*wt, 's', &
     acct_wall_time(qfall_acct)*tw, '%',                                &
    'Bound.conditions:',acct_cpu_time(bc_acct),                        's', &
     acct_cpu_time(bc_acct)*tc,    '%', acct_wall_time(bc_acct)*wt,    's', &
     acct_wall_time(bc_acct)*tw,  '%',                                  &
    'Message passing :',acct_cpu_time(mp_acct),                        's', &
     acct_cpu_time(mp_acct)*tc,    '%', acct_wall_time(mp_acct)*wt,    's', &
     acct_wall_time(mp_acct)*tw,  '%',                                  &
    'Miscellaneous   :',acct_cpu_time(misc_acct),                      's', &
     acct_cpu_time(misc_acct)*tc,  '%', acct_wall_time(misc_acct)*wt,  's', &
     acct_wall_time(misc_acct)*tw,'%'

    WRITE(6,'(/1x,a,2x,1p,e13.6,a,1x,0p,f7.2,a,1x,1p,e13.6,a,1x,0p,f7.2,a)')&
       'Entire model    :',total_cpu,'s',total_cpu/60.,' min',              &
                           total_wall*wt,'s',total_wall*wt/60.,' min'

    WRITE(6,'(/1x,a,2x,1p,e13.6,a,2x,0p,7x,4x,1p,e13.6,a,2x,0p,6x)')        &
       'Without Init/IO :',                                                 &
       total_cpu-acct_cpu_time(init_acct)-acct_cpu_time(output_acct),'s',   &
       (total_wall-acct_wall_time(init_acct)-acct_wall_time(output_acct))*wt,'s'

  END IF
  CALL flush(6)
  CALL mpbarrier

  RETURN
END SUBROUTINE acct_report_arps

SUBROUTINE init_walltime

  INTEGER :: lastval, rollover

  COMMON /walltime/ lastval, rollover

  lastval  = 0
  rollover = 0

  RETURN
END SUBROUTINE init_walltime

DOUBLE PRECISION FUNCTION f_walltime()
!
!-----------------------------------------------------------------------
!
!  Function to measure wall clock time.
!
!-----------------------------------------------------------------------
!
  INTEGER :: lastval, rollover
  COMMON /walltime/ lastval, rollover

  INTEGER :: wall_int
  INTEGER :: ticspersec,max_count

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL system_clock(wall_int,ticspersec,max_count)

  IF (wall_int < lastval .AND. lastval > max_count*0.01) THEN
    ! Assume that the clock has rolled over once.  The test for max_count*0.1
    ! is for Absoft which sometimes seems to have wall_int go backwards a
    ! little at times.
    WRITE (6,*) "F_WALLTIME: wall clock assumed to have rolled over once."
    rollover = rollover + 1
  END IF

  lastval = wall_int
  f_walltime =  dble(wall_int)/dble(ticspersec)                         &
             + (dble(max_count)/dble(ticspersec))*dble(rollover)

  RETURN
END FUNCTION f_walltime
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE FIX_SOIL_NSTYP               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE fix_soil_nstyp(nx,ny,nzsoil,nstypin,nstyp,                   &
                          tsoil,qsoil,wetcanp)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Adjust soil variables if nstyp differs between what is defined in
!  ARPS and what is in the data file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/08/11
!
!  MODIFICATION HISTORY:
!
!  ???
!  Soil variable update
!
!  12 June 2002 (Eric Kemp)
!  Bug fix for looping through tsoil,qsoil
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nstypin      Number of soil types in the data to be imported
!  nstyp        Number of soil types to use in the arps
!
!  INPUT/OUTPUT:
!
!  tsoil        Soil temperature (K)
!  qsoil        Soil moisture (m**3/m**3)
!  wetcanp      Canopy water amount
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx
  INTEGER :: ny
  INTEGER :: nzsoil
  INTEGER :: nstypin,nstyp

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyp) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyp) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,       0:nstyp) ! Canopy water amount

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Only one soil type, use average values.
!
!-----------------------------------------------------------------------

  IF ( nstypin == 1 ) THEN
    DO k = 1,nzsoil
      DO j = 1,ny-1
        DO i = 1,nx-1
          tsoil (i,j,k,1) = tsoil (i,j,k,0)
          qsoil (i,j,k,1) = qsoil (i,j,k,0)
        END DO
      END DO
    END DO
    DO j = 1,ny-1
      DO i = 1,nx-1
        wetcanp (i,j,1) = wetcanp (i,j,0)
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Use average values for soil types beyond that provided in the data.
!
!-----------------------------------------------------------------------

  IF ( nstypin < nstyp ) THEN
    DO is=nstypin+1,nstyp
      DO k = 1,nzsoil
        DO j = 1,ny-1
          DO i = 1,nx-1
            tsoil (i,j,k,is) = tsoil (i,j,k,0)
            qsoil (i,j,k,is) = qsoil (i,j,k,0)
          END DO
        END DO
      END DO
      DO j = 1,ny-1
        DO i = 1,nx-1
          wetcanp (i,j,is) = wetcanp (i,j,0)
        END DO
      END DO
    ENDDO
  END IF

!  GENE, WHAT IS THIS??????  6/03/02 DBW
!     tsoil (1:nx-1,1:ny-1,1:nzsoil-1,is) = tsoil (1:nx-1,1:ny-1, &
!            1:nzsoil-1,0)
!     qsoil (1:nx-1,1:ny-1,1:nzsoil-1,is) = qsoil (1:nx-1,1:ny-1, &
!            1:nzsoil-1,0)
!     wetcanp(1:nx-1,1:ny-1,is) = wetcanp(1:nx-1,1:ny-1,0)

  RETURN

END SUBROUTINE fix_soil_nstyp

!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE FIX_STYPFRCT_NSTYP             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE fix_stypfrct_nstyp(nx,ny,nstypin,nstyp,stypfrct)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Adjust stypfrct if the number of soil types read in differs from
!  the number specified for arps to use, nstyp.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/08/11
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nstypin      Number of soil types in the data file
!  nstyp        Number of soil types to use in the arps
!
!  INPUT/OUTPUT:
!
!  stypfrct     Soil type fraction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx
  INTEGER :: ny
  INTEGER :: nstypin,nstyp
  REAL    :: stypfrct(nx,ny,nstyp)

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,is
  REAL :: frctot

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

  INCLUDE "mp.inc"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Readjust stypfrct if only part of the soil type were read in.
!
!-----------------------------------------------------------------------
!
! 06/18/2002 Zuwen He
!
! Note: no exception handling for the nstypin == 0 case.
!
! There are also some other situations which should have been considered.
!

  IF (nstypin <= 0) THEN

    WRITE (6, '(1x,a)') ' ERROR: nstypin < 0, no way to fix it '
    WRITE (6, '(1x,a)') ' Program aborted in FIX_STYPFRCT_NSTYP.'
    CALL arpsstop('Stopped at calling fix_stypfrct_nstyp.',1)

  ELSE

!
! 06/18/2002 Zuwen
!
! otherwise, proceed ...
!
    IF (nstypin > nstyp) THEN

      IF (myproc == 0) WRITE(6,'(/,2(a,I2,a),/, 7x,2(a,I2),a,/)')       &
        ' INFO: read in "nstyp" from data file is ',nstypin,'.',        &
        ' The ARPS model runtime requires nstyp = ',nstyp,'.',          &
        'The extra data for "soiltyp" and "stypfrct" in section (:,:,', &
        nstyp+1,':',nstypin,') is discarded.'

      DO j=1,ny
        DO i=1,nx
          frctot = 0.0
          DO is=1,nstyp
            frctot = frctot + stypfrct(i,j,is)
          END DO
          IF (frctot == 0) THEN
            frctot = 1.
            stypfrct(i,j,1) = 1.
          ENDIF
          DO is=1,nstyp
            stypfrct(i,j,is) = stypfrct(i,j,is) / frctot
          END DO
        END DO
      END DO

    ELSE IF ( nstypin < nstyp ) THEN

!-----------------------------------------------------------------------
!
! Fill in any unset values.
!
!-----------------------------------------------------------------------

      IF (myproc == 0) WRITE(6,'(/,2(a,I2,a),/,10x,2(a,I2),a,/,10x,a,/)') &
        ' WARNING: read in "nstyp" from data file is ',nstypin,'.',       &
        ' The ARPS model runtime requires nstyp = ',nstyp,'.',            &
        'The 3D INTEGER array "soiltyp" will be zeros in section (:,:,',  &
        nstypin+1,':',nstyp,').',                                         &
        'You may encounter a floating exception in subroutine sfcphysics later.'

      stypfrct(1:nx-1,1:ny-1,nstypin+1:nstyp) = 0.

    END IF

!-----------------------------------------------------------------------
!
! Make sure the the fraction totals up to 1 (or 0)
!
!-----------------------------------------------------------------------

    DO j=1,ny
      DO i=1,nx
        frctot = 0.
        DO is = 1,nstyp
          frctot = frctot + stypfrct(i,j,is)
        END DO
        IF (frctot /= 0) THEN
          DO is = 1,nstyp
            stypfrct(i,j,is) = stypfrct(i,j,is) / frctot
          END DO
        END IF
      END DO
    END DO

  END IF !  nstyp <= 0

  RETURN

END SUBROUTINE fix_stypfrct_nstyp

!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE REMAP_SOIL_VARS              ######
!######                                                      ######
!######                    Developed by                      ######
!######            Weather Decision Technologies             ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE remap_soil_vars(nx,ny,nzsoil,nstypin,nstyp,                  &
                  tsoil_in,qsoil_in,wetcanp_in,soiltyp_in,              &
                  tsfcin,tsoilin,wsfcin,wdpin,qsoilin,wcanpin,          &
                  intver,                                               &
                  tsoil,qsoil,wetcanp,soiltyp)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Rearrange soil variables to be consistent with any difference between
!  the soil types to be used and the ones in the data set.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/10/27
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx
!  ny
!  nzsoil
!  nstypin      Number of soil types in the data file
!  nstyp        Number of soil types to use in the arps
!  tsoilin
!  qsoilin
!  wcanpin
!  tsoil_in     Soil temperature in data set (K)
!  qsoil_in     Soil moisture in data set
!  wetcanp_in   Canopy water amount in data set
!  soiltyp_in   Soil type in data set
!
!  OUTPUT:
!
!  tsoil        Soil temperature (K)
!  qsoil        Soil moisture (m**3/m**3)
!  wetcanp      Canopy water amount
!  soiltyp      Soil type
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx, ny, nzsoil
  INTEGER :: nstypin,nstyp
  INTEGER :: tsfcin,tsoilin,wsfcin,wdpin,qsoilin,wcanpin

  REAL    :: tsoil_in (nx,ny,nzsoil,0:nstypin) ! Soil temperature (K)
  REAL    :: qsoil_in (nx,ny,nzsoil,0:nstypin) ! Soil moisture (m**3/m**3)
  REAL    :: wetcanp_in(nx,ny,0:nstypin)       ! Canopy water amount
  INTEGER :: soiltyp_in(nx,ny,nstypin)         ! Soil type in model domain

  REAL :: tsoil (nx,ny,nzsoil,0:nstyp) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyp) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyp)       ! Canopy water amount
  INTEGER :: soiltyp(nx,ny,nstyp)      ! Soil type in model domain

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,is,isin
  INTEGER :: ifound_styp
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
!  CHARACTER (LEN=40) :: fmtver,fmtver410,fmtver500
!  PARAMETER (fmtver410='* 004.10 GrADS Soilvar Data',intver410=410)
!  PARAMETER (fmtver500='* 005.00 GrADS Soilvar Data',intver500=500)

  INTEGER, PARAMETER  :: intver410 = 410, intver500 = 500

  INTEGER, INTENT(IN) :: intver

!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Match the soil type in the old data set.  If not found then set to
! average soil properties (except water).
!
!-----------------------------------------------------------------------

  IF (intver == intver410) THEN

    DO j=1,ny
      DO i=1,nx
        DO is=1,nstyp

        ifound_styp = 0
        DO isin=1,nstypin
          IF (soiltyp(i,j,is) == soiltyp_in(i,j,isin)) THEN
            IF (tsfcin == 1)  tsoil(i,j,1,is) = tsoil_in(i,j,1,isin)
            IF (tsoilin == 1) tsoil(i,j,2,is) = tsoil_in(i,j,2,isin)
            IF (wsfcin == 1)  qsoil(i,j,1,is) = qsoil_in(i,j,1,isin)
            IF (wdpin == 1)   qsoil(i,j,2,is) = qsoil_in(i,j,2,isin)
            IF (wcanpin == 1) wetcanp(i,j,is) = wetcanp_in(i,j,isin)
            ifound_styp = 1
            EXIT
          END IF
        END DO
        IF (ifound_styp == 0) THEN
          IF (tsfcin == 1)  tsoil(i,j,1,is) = tsoil_in(i,j,1,0)
          IF (tsoilin == 1) tsoil(i,j,2,is) = tsoil_in(i,j,2,0)
          IF (wsfcin == 1)  qsoil(i,j,1,is) = qsoil_in(i,j,1,0)
          IF (wdpin == 1)   qsoil(i,j,2,is) = qsoil_in(i,j,2,0)
          IF (wcanpin == 1) wetcanp(i,j,is) = wetcanp_in(i,j,0)
          IF (tsfcin == 1 .and. soiltyp(i,j,is) == 13) THEN
            tsoil(i,j,1,is) = tsoil_in(i,j,1,0)
          END IF
        END IF

        END DO
      END DO
    END DO

  ELSE IF (intver >= intver500) THEN

    DO k=1,nzsoil
      DO j=1,ny
        DO i=1,nx
          DO is=1,nstyp

          ifound_styp = 0
          DO isin=1,nstypin
            IF (soiltyp(i,j,is) == soiltyp_in(i,j,isin)) THEN
              IF (tsoilin == 1) tsoil(i,j,k,is) = tsoil_in(i,j,k,isin)
              IF (qsoilin == 1) qsoil(i,j,k,is) = qsoil_in(i,j,k,isin)
!
! 07/08/2002 Zuwen He
!
! wetcanp only need to be assigned once.
!
              IF (k == 1 .AND. wcanpin == 1) wetcanp(i,j,is)=wetcanp_in(i,j,isin)
              ifound_styp = 1
              EXIT
            END IF
          END DO
          IF (ifound_styp == 0) THEN
            IF (tsoilin == 1) tsoil(i,j,k,is) = tsoil_in(i,j,k,0)
            IF (qsoilin == 1) qsoil(i,j,k,is) = qsoil_in(i,j,k,0)
            IF (k == 1 .AND. wcanpin == 1) wetcanp(i,j,is) = wetcanp_in(i,j,0)
            IF (tsoilin == 1 .and. soiltyp(i,j,is) == 13) THEN
              tsoil(i,j,k,is) = tsoil_in(i,j,k,0)
            END IF
          END IF

          END DO
        END DO
      END DO
    END DO

  ELSE

!-----------------------------------------------------------------
!
! 07/08/2002  Zuwen He
!
! although it wont happen, ...
!
!----------------------------------------------------------------

    write (6,*) "Incorrect vertion in remap_soil_vars"
    CALL arpsstop("arpsstop called from remap_soil_vars",1)

  END IF ! intver

END SUBROUTINE remap_soil_vars


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LD2DGRID                   ######
!######                                                      ######
!######                    Developed by                      ######
!######            Weather Decision Technologies             ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE ld2dgrid(dxin,dyin,ctrlatin,ctrlonin,  &
           mprojin,trlat1in,trlat2in,trlonin,sclfctin)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Load 2-D grid parameters into the ARPS common blocks (useful for
!  C programs using ARPS subroutines).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2001/08/03
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL :: dxin,dyin,ctrlatin,ctrlonin,trlat1in,trlat2in,trlonin,sclfctin
  INTEGER :: mprojin

!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  dx = dxin
  dy = dyin
  ctrlat = ctrlatin
  ctrlon = ctrlonin
  mapproj = mprojin
  trulat1 = trlat1in
  trulat2 = trlat2in
  trulon = trlonin
  sclfct = sclfctin

  RETURN
END SUBROUTINE ld2dgrid

!
! Convert a character string to upper case
!
FUNCTION upcase(string) RESULT(upper)

  IMPLICIT NONE

  INTEGER, PARAMETER :: lenstr = 40

  CHARACTER(LEN=lenstr), INTENT(IN) :: string
  CHARACTER(LEN=lenstr)             :: upper

  INTEGER :: j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  upper = ' '
  DO j = 1,lenstr
    IF(string(j:j) >= "a" .AND. string(j:j) <= "z") THEN
      upper(j:j) = ACHAR(IACHAR(string(j:j)) - 32)
    ELSE
      upper(j:j) = string(j:j)
    END IF
  END DO

END FUNCTION upcase

!WDT testing:

!SUBROUTINE test_dump(data,filename,nx,ny,nz,stagdim,incfake)
!
!      ! write out a unformatted dump to debug MP version
!
!      implicit none
!
!      include 'mp.inc'
!
!      integer nx,ny,nz
!      integer stagdim,incfake
!      real data(nx,ny,nz)
!      character*(*) filename
!      character*80 tempname
!
!      integer count
!      common /save_count/ count
!
!      !return
!      count = count + 1
!
!      write(tempname, '(a,a,i4.4)') trim(filename),'_c',count
!
!      if (mp_opt > 0) then
!        write(tempname, '(a,a,2i2.2)')  &
!           trim(adjustl(tempname)),'_',loc_x,loc_y
!      endif
!
!      open (11,file=tempname,form='unformatted')
!      write (11) nx,ny,nz,stagdim,incfake
!      write (11) data
!      close (11)
!
!      CALL test_check_in(tempname)
!
!      RETURN
!END
!
!SUBROUTINE test_dump2(data,filename,nx,ny,nz,ib,ie,jb,je,kb,ke)
!
!      ! write out a unformatted dump to debug MP version
!
!      implicit none
!
!      include 'mp.inc'
!
!      integer nx,ny,nz
!      integer ib,ie,jb,je,kb,ke   ! compare for i=ib,ie ...
!      real data(nx,ny,nz)
!      character*(*) filename
!      character*80 tempname
!
!      integer count
!      common /save_count/ count
!
!      !return
!      count = count + 1
!
!      write(tempname, '(a,a,i4.4)') trim(filename),'_c',count
!
!      if (mp_opt > 0) then
!        write(tempname, '(a,a,2i2.2)')  &
!           trim(adjustl(tempname)),'_',loc_x,loc_y
!      endif
!
!      open (11,file=tempname,form='unformatted')
!      write (11) nx,ny,nz,-1,-1
!      write (11) ib,ie,jb,je,kb,ke
!      write (11) data
!      close (11)
!
!      RETURN
!END
!
!subroutine test_check_in(comment)
!      implicit none
!      include 'mp.inc'
!      character*(*) comment
!
!      !return
!      CALL flush(6)
!      CALL flush(0)
!      write (*,*) "XXX checking in:",myproc," ",trim(comment)
!      CALL flush(6)
!      CALL mpbarrier
!
!      return
!end
!!!
