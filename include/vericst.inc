!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                    VERICST.INC                       ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     Center for Analysis and Prediction of Storms     ######
!     ######                University of Oklahoma                ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
!
!#######################################################################
!
!     PURPOSE:
!
!     Include file 'vericst.inc' for ARPS verification software.
!
!     This file defines the arrays to be used in the verification
!     code for the ARPS model.
!
!#######################################################################
!
!     AUTHOR: John Mewes
!     09/15/97
!
!     MODIFICATION HISTORY:
!
!     Eric Kemp, May 2000
!     Added precipitation station arrays.
!
!     Eric Kemp, June 2000
!     Added MOS station arrays.
!
!#######################################################################
!
!     Sounding verification station data variables
!
!#######################################################################
!
      integer, parameter :: sndmax=2000 ! Maximum number of stations in sndstns.tbl
      integer :: sndstn                 ! Actual number of stations in the domain

      character(4) :: sndstid(sndmax) ! Sounding stn's letter identifier

      real :: sndlat(sndmax)         ! Station latitude
      real :: sndlon(sndmax)         ! Station longitude
      real :: sndelev_act(sndmax)    ! Station elevation (actual)

      real :: sndxpt(sndmax),sndypt(sndmax)    ! Stn locations in x,y coords
      integer :: sndipt(sndmax),sndjpt(sndmax) ! Stn locations in i,j coords

      common /sndver1/ sndstn,sndlat,sndlon,sndelev_act,&
&sndxpt,sndypt,sndipt,sndjpt
      common /sndver2/ sndstid

!#######################################################################
!
!     Profiler verification station data variables
!
!#######################################################################

      integer, parameter :: promax=2000  ! Maximum number of stations in profstns.tbl
      integer :: prostn     ! Actual number of stations in the domain

      character(3) :: prostid(promax) ! Profiler stn's letter identifier

      real :: prolat(promax)         ! Station latitude
      real :: prolon(promax)         ! Station longitude
      real :: proelev_act(promax)    ! Station elevation (actual)

      real :: proxpt(promax),proypt(promax)    ! Stn locations in x,y coords
      integer :: proipt(promax),projpt(promax) ! Stn locations in i,j coords

      common /prover1/ prostn,prolat,prolon,proelev_act,&
&proxpt,proypt,proipt,projpt
      common /prover2/ prostid

!#######################################################################
!
!     Surface verification station data variables
!
!#######################################################################

      integer, parameter :: sfcmax=5000  ! Maximum number of stations in msstns.tbl
      integer :: sfcstn     ! Actual number of stations in the domain
      integer :: sfcstn_master ! Master list

      character(4) ::  sfcstid(sfcmax) ! Surface stn's letter identifier
      character(4) ::  sfcstid_master(sfcmax) ! Master list
      character(4) ::  sfcuse(sfcmax)  ! Subset we want to use

      real :: sfclat(sfcmax)         ! Station latitude
      real :: sfclon(sfcmax)         ! Station longitude
      real :: sfcelev_act(sfcmax)    ! Station elevevation (actual)

      real :: sfcxpt(sfcmax),sfcypt(sfcmax)    ! Stn locations in x,y coords
      integer :: sfcipt(sfcmax),sfcjpt(sfcmax) ! Stn locations in i,j coords

      integer :: nsfcuse

      logical :: readstns

      common /sfcver1/ sfcstn,sfcstn_master,sfclat,sfclon,sfcelev_act,&
&sfcxpt,sfcypt,sfcipt,sfcjpt,readstns
      common /sfcver2/ sfcstid, sfcstid_master
      common /sfcver3/ nsfcuse, sfcuse


!#######################################################################
!
!     Precipitation verification station data variables
!
!#######################################################################

      integer, parameter :: precmax=2000  ! Maximum number of stations in msstns.tbl
      integer :: precstn     ! Actual number of stations in the domain

      character(16) ::  precstid(precmax) ! Surface stn's letter identifier

      real :: preclat(precmax)         ! Station latitude
      real :: preclon(precmax)         ! Station longitude
      real :: precelev_act(precmax)    ! Station elevevation (actual)

      real :: precxpt(precmax),precypt(precmax)    ! Stn locations in x,y
                                                   ! coords
      integer :: precipt(precmax),precjpt(precmax) ! Stn locations in i,j
                                                   ! coords

      common /precver1/ precstn,preclat,preclon,precelev_act,&
&precxpt,precypt,precipt,precjpt
      common /precver2/ precstid

!#######################################################################
!
!     MOS verification station data variables
!
!#######################################################################

      integer, parameter :: mosmax=2000   ! Maximum number of stations in mosstns.tbl
      integer :: mosstn     ! Actual number of stations in the domain

      character(4) :: mosstid(mosmax) ! Sounding stn's letter identifier

      real :: moslat(mosmax)         ! Station latitude
      real :: moslon(mosmax)         ! Station longitude
      real :: moselev_act(mosmax)    ! Station elevation (actual)

      real :: mosxpt(mosmax),mosypt(mosmax)    ! Stn locations in x,y coords
      integer :: mosipt(mosmax),mosjpt(mosmax) ! Stn locations in i,j coords

      common /mosver1/ mosstn,moslat,moslon,moselev_act,&
&mosxpt,mosypt,mosipt,mosjpt
      common /sndver2/ mosstid

