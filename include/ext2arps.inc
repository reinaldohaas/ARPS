!
!-----------------------------------------------------------------------
!
!  NAMELIST parameters for EXT2ARPS (In arps.input)
!
!-----------------------------------------------------------------------
!
  INTEGER :: iorder              ! Order of polynomial used for interpolation
                                 ! = 1  Linear
                                 ! = 2  Quadratic
                                 ! = 3  Cubic
  INTEGER :: intropt             ! Option indicating to interpolate
                                 ! perturbation or total variables:
                                 ! = 1  Interpolate perturbation variables
                                 !      and add to base sounding (default);
                                 ! = 2  Interploate total variables (except
                                 !      pressure).
  INTEGER :: nsmooth             ! Number of 27-pt (3-D 3 pt)
                                 ! smoothing passes after interpolation.
                                 ! 1 or 2 recommended.
  INTEGER :: ext_lbc             ! Option to apply lateral boundary conditions
                                 ! to the winds.
                                 ! = 0  no boundary contitions applied;
                                 ! = 1  apply zero-gradient boundary contitions
                                 !      (default).
  INTEGER :: ext_vbc             ! Option to apply vertical boundary conditions
                                 ! to w.
                                 ! = 0  no boundary contitions applied;
                                 ! = 1  apply boudary contitions specified
                                 !      by tbc & bbc to w (default).
  INTEGER :: exttrnopt           ! Terrain option for output grid:
                                 ! = 0  use terrain specified by ternopt
                                 !      above (default);
                                 ! = 1  interpolate terrain from original grid.

  INTEGER :: extsfcopt           ! Option for near surface interpolation:
                                 ! = 0  use pressure level fields (default)
                                 ! = 1  use near surface fields available
                                 !      (simple interpolation)

  INTEGER :: extntmrg            ! Number of zones to merge for exttrnopt=2


  CHARACTER (LEN=256) :: dir_extd ! directory of external data

  CHARACTER (LEN=256) :: extdname ! Prefix string of external file name.
                                  ! Used ONLY for ARPS data input
                                  ! The file name should look like:
                                  !      extdname.time_string

  INTEGER :: extdopt             ! External data source option
                                 ! = 0, ARPS, default
                                 ! = 1, NMC RUC (Hybrid) data
                                 ! = 2, NMC ETA data
                                 ! >= 3, User defined
                                 !       (user should provide programs)

  INTEGER :: extdfmt             ! External data format indicator. For
                                 ! extdopt=0, see description of hdmpfmt

  INTEGER :: nextdfil            ! number of external data files to process

  INTEGER :: grdbasopt           ! grdbas history file option
                                 ! = -1, never write files (LBC run).  Output
                                 !       from a "t0" run is sufficient.
                                 ! =  0, for LBC runs, only the first "extm"
                                 !       file produces output.  This is a big
                                 !       I/O saver if the files aren''t used.
                                 ! =  1, always write the "grdbas" files.

  CHARACTER (LEN=29) :: extdtime(50)   ! external data times
                                       ! format mm-dd-yyyy:hh:mm:ss/HHH:MM:SS
                                       !        _initial time_____ /forecast

  INTEGER :: i2dfmt              ! 2D field output option
                                 ! = 0,  no 2D field output
                                 ! = 1,  native binary 2D output
                                 ! = 3,  HDF4 2D output
                                 ! = 7,  NetCDF 2D output

  CHARACTER (LEN=256) :: outheader     ! header of 2D output file name

  INTEGER :: iboxs, iboxe, jboxs, jboxe

  NAMELIST /extdfile/ extdopt,extdfmt,dir_extd,extdname,nextdfil,extdtime,  &
                      iorder,nsmooth,ext_lbc,ext_vbc,intropt,exttrnopt, &
                      extntmrg,extsfcopt,grdbasopt,i2dfmt,outheader,    &
                      iboxs, iboxe, jboxs, jboxe
