!-----------------------------------------------------------------------
!  Include file 'arpssfc.inc' for ARPS arpssfc.f program
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  1/26/94
!
!  MODIFICATIONS:
!
!  04/15/98 (Dan Weber)
!  Added gvnx,gvny, and gvfrcopt for use with the NESDIS green
!  vegetation data set.
!
!  2000/08/23 (Gene Bassett)
!  Moved the data set parameter statements to namelist input file.
!
!  2001/07/26 (Gene Bassett)
!  Eliminated variables colst,rowst,colvt,rowvt,colnv,&rownv.
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: gvnx = 2500  ! Number of NESDIS green veg. data points (e-w)
  INTEGER, PARAMETER :: gvny =  904  ! Number of NESDIS green veg. data points (n-s)

  INTEGER, PARAMETER :: nsoiltyp=13  ! number of soil types
  INTEGER, PARAMETER :: nvegtyp =14  ! number of vegatation types
!
!-----------------------------------------------------------------------
!
!  Data flags read in from input namelist &outflag
!
!      stypout: 1/0 - yes/no, Output styp
!      vtypout: 1/0 - yes/no, Output vtyp
!       laiout: 1/0 - yes/no, Output lai
!      rfnsout: 1/0 - yes/no, Output roufns
!
!-----------------------------------------------------------------------
!
  INTEGER :: stypout,vtypout,laiout,rfnsout,vegout,ndviout

  COMMON /comflag/ stypout,vtypout,laiout,rfnsout,vegout,ndviout
