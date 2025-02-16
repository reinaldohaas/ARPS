!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########        ENSEMBLE        ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################
!
!
!
!======================================================================-
!             &&&&    G E N E R A L    I N F O R M A T I O N    &&&&
!======================================================================-
!
!  AUTHOR:  Henry Neeman, Univ. of Oklahoma, Project SAMEX
!           (hneeman@tornado.caps.ou.edu)
!
!  HISTORY:
!        1998/03/18  First written [HJN].
!        1999/11/19  Output grid parameters moved to input file [Gene Bassett].
!
      ! Maximum size of 2-D output grid.
  INTEGER :: max2d_out
  PARAMETER (max2d_out=256*256)  ! MAX2D_OUT needs to be >= qni*qnj
  INTEGER :: qni, qnj               ! output grid dimensions, set in input file
  COMMON /outconfigi/ qni, qnj
  CHARACTER (LEN=256) :: enstag, instit, model
  COMMON /outconfigc/ enstag, instit, model
  REAL :: qswcorn, qswcorw, qnecorn, qnecorw
  COMMON /outconfigr/ qswcorn, qswcorw, qnecorn, qnecorw

  CHARACTER (LEN=*) :: resun
  PARAMETER (resun='km')
  CHARACTER (LEN=*) :: mslptag, mslpunits
  PARAMETER (mslptag='MSLP', mslpunits='pascals')
  CHARACTER (LEN=*) :: hgt250tag, hgt250units
  PARAMETER (hgt250tag='HGT250', hgt250units='m')
  CHARACTER (LEN=*) :: vort250tag, vort250units
  PARAMETER (vort250tag='ABSVORT250', vort250units='1/s')
  CHARACTER (LEN=*) :: uwind250tag, uwind250units
  PARAMETER (uwind250tag='UWIND250', uwind250units='m/s')
  CHARACTER (LEN=*) :: vwind250tag, vwind250units
  PARAMETER (vwind250tag='VWIND250', vwind250units='m/s')
  CHARACTER (LEN=*) :: hgt500tag, hgt500units
  PARAMETER (hgt500tag='HGT500', hgt500units='m')
  CHARACTER (LEN=*) :: vort500tag, vort500units
  PARAMETER (vort500tag='ABSVORT500', vort500units='1/s')
  CHARACTER (LEN=*) :: uwind500tag, uwind500units
  PARAMETER (uwind500tag='UWIND500', uwind500units='m/s')
  CHARACTER (LEN=*) :: vwind500tag, vwind500units
  PARAMETER (vwind500tag='VWIND500', vwind500units='m/s')
  CHARACTER (LEN=*) :: hgt850tag, hgt850units
  PARAMETER (hgt850tag='HGT850', hgt850units='m')
  CHARACTER (LEN=*) :: vort850tag, vort850units
  PARAMETER (vort850tag='ABSVORT850', vort850units='1/s')
  CHARACTER (LEN=*) :: uwind850tag, uwind850units
  PARAMETER (uwind850tag='UWIND850', uwind850units='m/s')
  CHARACTER (LEN=*) :: vwind850tag, vwind850units
  PARAMETER (vwind850tag='VWIND850', vwind850units='m/s')
  CHARACTER (LEN=*) :: temp850tag, temp850units
  PARAMETER (temp850tag='TEMP850', temp850units='degC')
  CHARACTER (LEN=*) :: thcktag, thckunits
  PARAMETER (thcktag='THCK', thckunits='m')
  CHARACTER (LEN=*) :: rh700tag, rh700units
  PARAMETER (rh700tag='RH700', rh700units='%')
  CHARACTER (LEN=*) :: omega700tag, omega700units
  PARAMETER (omega700tag='OMEGA700', omega700units='-microb/s')
  CHARACTER (LEN=*) :: temp2mtag, temp2munits
  PARAMETER (temp2mtag='TEMP2M', temp2munits='degC')
  CHARACTER (LEN=*) :: dewp2mtag, dewp2munits
  PARAMETER (dewp2mtag='DEWP2M', dewp2munits='degC')
  CHARACTER (LEN=*) :: accppttag, accpptunits
  PARAMETER (accppttag='ACCPPT', accpptunits='mm')
  CHARACTER (LEN=*) :: convppttag, convpptunits
  PARAMETER (convppttag='CONVPPT', convpptunits='mm')
  CHARACTER (LEN=*) :: srehtag, srehunits
  PARAMETER (srehtag='SREH', srehunits='m/s^2')
  CHARACTER (LEN=*) :: litag, liunits
  PARAMETER (litag='LI', liunits='degC')
  CHARACTER (LEN=*) :: capetag, capeunits
  PARAMETER (capetag='CAPE', capeunits='J/kg')
  CHARACTER (LEN=*) :: brntag, brnunits
  PARAMETER (brntag='BRN', brnunits='[unitless]')
  CHARACTER (LEN=*) :: cintag, cinunits
  PARAMETER (cintag='CIN', cinunits='J/kg')
  CHARACTER (LEN=*) :: llmctag, llmcunits
  PARAMETER (llmctag='LLMC', llmcunits='g/kg/s')
  CHARACTER (LEN=*) :: pwtag, pwunits
  PARAMETER (pwtag='PW', pwunits='g/m^2')
  CHARACTER (LEN=*) :: crftag, crfunits
  PARAMETER (crftag='COMPREFL', crfunits='dbz')
  CHARACTER (LEN=*) :: qstag
  PARAMETER (qstag='CENTER')

