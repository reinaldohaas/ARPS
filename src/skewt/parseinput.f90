
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########       PARSEINPUT       ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################


SUBROUTINE ParseInput (                                                 &
      nsmax         , nsnd          , sndfile       , has_wind      ,   &
      print_info    , verbose       , logfiles      , do_indices    ,   &
      color_file    , outfile       ,                                   &
      do_hodo       , hodo_denmwind , jhunits       , helcontrs     ,   &
      plot_cb_parcel,plot_sfc_parcel, plot_tv       ,                   &
      modify_snd    , weightq       , write_modsnd  ,                   &
      tv_wl         , cape_use_t    , cape_use_irrev,                   &
      tmin          , tmax          , pmin          , pmax          )

!======================================================================-
!             &&&&    G E N E R A L    I N F O R M A T I O N    &&&&
!======================================================================-
!
!  PURPOSE:  This subroutine parses command-line arguments for the Skew-T
!                program.
!
!  AUTHOR:  Richard Carpenter, Univ. of Oklahoma
!            (rcarpenter@ou.edu)
!
!  HISTORY:
!        09/13/93  Created from u.skewt.f
!        1995/01/24  -tvnowl now also plots Tv of environment
!
!  INPUT ARGUMENTS: <none>
!
!  IN/OUT ARGUMENTS: <none>
!
!  OUTPUT ARGUMENTS:  [.] means logical, [$] means CHARACTER, (,) means array.
!    Variable                Cmd-line flag        Meaning
!    --------                -------------        -------
!    verbose[.]                -v                Be verbose
!    do_hodo[.]                -hodo                Also plot hodo
!    helcontrs[.]        -helcontrs        Plot helicity contours on hodograph.
!    plot_sfc_parcel[.]        -sfc                Plot parcel ascending from sfc
!    plot_cb_parcel[.]        -cb                Plot parcel ascending from
!                                        specified cloud base
!    plot_tv[.]                -tv,-tvnowl        Also plot Tv
!    tv_wl[.]                -tv,-tvnowl        Use water loading when plotting Tv.
!    has_wind[.]        -haswind        Sounding file also has wind.
!    modify_snd[.]        -mod                CALL Keith Brewster's routine
!                                        to modify sounding
!    weightq[.]                -weightq        In KB's routine, modify Qv
!                                        using a weighting function
!    write_modsnd        -writemod        Write the modified sounding to file.mod
!    cape_use_t[.]        -capet                Compute CAPE using T not Tv
!    cape_use_irrev[.]        -capeirrev        Compute CAPE for pseudoadiab,
!                                        not reversible, ascent
!    print_info[.]        -info                Forces printng (on plot) of certain info
!    nsnd                <n/a>                Number of sounding files
!    sndfile[$]                [file]                Main sounding file
!    sndfile2[$]()        [file(s)]        Optional additional snd files
!        <n/a>                -jerry                -writemod -sfc -weightq -capet
!                                        -capeirrev
!    tmin,tmax,pmin,pmax: -coords tmin,tmax,pmin,pmax
!                                        Use to window in on a region. Units are
!                                        K, mb.
!
!  I/O:        <none>
!
!  SPECIAL REQUIREMENTS:
!        IARGC, GETARGS
!  
!  OTHER INFORMATION:
!        A parallel version, for ARPSTools, reads a NAMELIST instead.
!
!======================================================================-
!                    %%%%    D E C L A R A T I O N S    %%%%
!======================================================================-

  IMPLICIT NONE

! - - - - - - - - -  Include files - - - - - - - - - - - - - - - - - - -
! - - - - - - - - -  Constant declarations - - - - - - - - - - - - - - -
! - - - - - - - - -  Argument declarations - - - - - - - - - - - - - - -

  INTEGER, INTENT(IN)  :: nsmax
  REAL,    INTENT(OUT) :: tmin, tmax, pmin, pmax
  INTEGER, INTENT(OUT) :: nsnd, jhunits
  LOGICAL, INTENT(OUT) ::                                               &
      do_hodo       , verbose       , has_wind      , logfiles      ,   &
      plot_cb_parcel,plot_sfc_parcel, plot_tv       , tv_wl         ,   &
      print_info    , do_indices    ,                                   &
      modify_snd    , weightq       , write_modsnd  , helcontrs     ,   &
      hodo_denmwind , cape_use_t    , cape_use_irrev

  CHARACTER (LEN=*), INTENT(OUT) :: color_file, sndfile(nsmax), outfile

! - - - - - - - - -  Global/External declarations  - - - - - - - - - - - 

  INTEGER, EXTERNAL :: IARGC

! - - - - - - - - -  Local declarations  - - - - - - - - - - - - - - - - 

  INTEGER                i             , nargs
  REAL        save          
  CHARACTER*255        string
  CHARACTER*79        usagestr(36)
  DATA        usagestr / &
' ',        &
'   This program produces a Skew-T diagram from a raw or model sounding.', &
' ',        &
'Usage: skewt [options] sounding-file[s]',                              &
' ',        &
' Options:',                                                            &
'   -v:         Verbose',                                               &
'   -hodo:      Plot Hodograph',                                        &
'   -helcontrs: Plot helicity contours on hodograph',                   &
'   -denmw:     Compute & plot density wgtd mean wind on hodo',         &
'   -sfc:       Plot T of parcel lifted from surface ',                 &
'   -indices:   Compute complete set of stability indices (implies -sfc)', &
'   -cb:        Plot T of parcel lifted from (specified) cloud base',   &
'   -tvwl:      Plot Tv (w/ water loading) instead of T for parcels;',  &
'               also plot Tv of sounding',                              &
'   -tv:        Plot Tv (w/o water loading) instead of T for parcels;', &
'               also plot Tv of sounding',                              &
'   -haswind:   A raw sounding file has the format p,T,Td,u,v',         &
'   -mod:       Modify sounding using Keith Brewster''s routine',       &
'   -weightq:   Weight mixing ratio when modifying sounding',           &
'   -writemod:  Write the modified sounding to file.mod',               &
'   -capet:     Use T instead of Tv when computing CAPE',               &
'   -caperev:   Use reversible adiabat',                                &
'   -capeirrev: Use pseudo- rather than reversible adiabat',            &
'   -info:      Write parcel info on graph even if not plotting parcel',&
'   -2:         Plot additional info on a second graph',                &
'   -coords tmin tmax pmin pmax: Zoom in on a region Units: (K, mb)',   &
'   -jerry:     Apply options -sfc, -weightq, -capet, -capeirrev',      &
'   -standard:  Apply options -hodo, -sfc, -hunits 1',                  &
'   -logfiles:  Output soundings and indices',                          &
'   -hunits:    0=km, 1=kft',                                           &
'   -colorfile file: File containing colors for graphical elements',    &
'   -o file    output file name',                                       &
' ' ,                                                                   &
' ARPS format sounding: First character of file must be "&".' ,         &
' ' /

!======================================================================-
!                 @@@@    E X E C U T A B L E    C O D E    @@@@
!======================================================================-

! First get environmental variables

  CALL GETENV ('SKEWT_COLORFILE', string)
  IF (LEN_TRIM(string) > 0) THEN
    color_file = TRIM(ADJUSTL(string))
    PRINT *, 'Using environmental variable SKEWT_COLORFILE=', TRIM(color_file)
  END IF

  nsnd        = 0

  nargs = IARGC()
  IF (nargs.EQ.0) THEN
    PRINT '(1x,a)', usagestr
    STOP
  END IF

  i = 0
  DO WHILE (i.LT.nargs)
    i        = i + 1
    CALL GETARG (i,string)

  ! Read option

    IF (string(1:1).EQ.'-') THEN

      IF (verbose) PRINT *, 'Option: ', TRIM(string)

      IF (string.EQ.'-coords') THEN
        i = i + 1
        CALL GETARG (i,string)
        READ (string,*) tmin
        i = i + 1
        CALL GETARG (i,string)
        READ (string,*) tmax
        i = i + 1
        CALL GETARG (i,string)
        READ (string,*) pmin
        i = i + 1
        CALL GETARG (i,string)
        READ (string,*) pmax
        IF (tmin.GT.tmax) THEN
          save        = tmax
          tmax        = tmin
          tmin        = save
        END IF
        IF (pmin.GT.pmax) THEN
          save        = pmax
          pmax        = pmin
          pmin        = save
        END IF
        PRINT '(1x,a,4f7.1)','Changing default coords to (T,p): ', &
           tmin,tmax,pmin,pmax

      ELSE IF (string(:6).EQ.'-color') THEN
          i = i + 1
          CALL GETARG (i,color_file)

      ELSE IF (string(:2).EQ.'-o') THEN
          i = i + 1
          CALL GETARG (i,outfile)

      ELSE IF (string.EQ.'-logfiles') THEN
          logfiles = .TRUE.

      ELSE IF (string.EQ.'-caperev') THEN
          cape_use_irrev = .FALSE.

      ELSE IF (string.EQ.'-capeirrev') THEN
          cape_use_irrev = .TRUE.

      ELSE IF (string.EQ.'-capet') THEN
          cape_use_t = .TRUE.

      ELSE IF (string.EQ.'-cb') THEN
          plot_cb_parcel = .TRUE.

      ELSE IF (string.EQ.'-denmw') THEN
          hodo_denmwind = .TRUE.

      ELSE IF (string.EQ.'-haswind') THEN
          has_wind = .TRUE.

      ELSE IF (string.EQ.'-helcontrs') THEN
          helcontrs = .TRUE.

      ELSE IF (string.EQ.'-hodo') THEN
          do_hodo = .TRUE.

      ELSE IF (string.EQ.'-hunits') THEN
          i = i + 1
          CALL GETARG (i,string)
          READ (string,*) jhunits

      ELSE IF (string.EQ.'-info') THEN
          print_info = .TRUE.
          plot_sfc_parcel = .TRUE.

      ELSE IF (string(:4).EQ.'-ind') THEN
          do_indices = .TRUE.

      ELSE IF (string.EQ.'-writemod') THEN
          write_modsnd        = .TRUE.

      ELSE IF (string.EQ.'-jerry') THEN
          write_modsnd        = .TRUE.
          plot_sfc_parcel = .TRUE.
          weightq = .TRUE.
          cape_use_t = .TRUE.
          cape_use_irrev = .TRUE.

      ELSE IF (string.EQ.'-std' .OR. string.EQ.'-standard') THEN
          jhunits = 1
          plot_sfc_parcel = .TRUE.
          do_hodo = .TRUE.
          do_indices = .TRUE.

      ELSE IF (string.EQ.'-mod') THEN
          modify_snd = .TRUE.

      ELSE IF (string.EQ.'-sfc') THEN
          plot_sfc_parcel = .TRUE.

      ELSE IF (string.EQ.'-tvwl') THEN
          plot_tv = .TRUE.
          tv_wl        = .TRUE.

      ELSE IF (string.EQ.'-tv') THEN
          plot_tv = .TRUE.
          tv_wl        = .FALSE.

      ELSE IF (string.EQ.'-v') THEN
          verbose = .TRUE.

      ELSE IF (string.EQ.'-weightq') THEN
          weightq = .TRUE.

      ELSE
          PRINT *, 'Unknown option: ', string
          PRINT '(1x,a)', usagestr
          STOP
      END IF

    ! READ file name

    ELSE
      nsnd        = nsnd + 1
      sndfile(nsnd)        = string
      PRINT *, 'Sounding file: ', TRIM(sndfile(nsnd))
    END IF

  END DO

END

!
! Added arpssighandle, used only on IBM p690 platform
!
! To avoid having to link with libarps.a
!
SUBROUTINE arpssighandle()

  IMPLICIT NONE

  WRITE(6,*) "ERROR:  caught a signal."
  WRITE(0,*) "arpssighandle called due to a float signal."
  STOP

  RETURN
END SUBROUTINE arpssighandle

