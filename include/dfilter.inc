
  integer, parameter ::  df_maxnstps = 100
  integer, parameter ::  numdfvar3d=8, numdfvar_soil=4,            &
                         numdfvar_soil2d=1, numdfvar_soil2ds=1

  INTEGER :: numdfvar_atm3d

  integer :: dfilter_opt
  real    :: df_tstart,df_tinv
  integer :: df_nstps
  real    :: df_wght(df_maxnstps)

  common /df_01/ dfilter_opt,df_tstart,df_tinv,df_nstps,df_wght

