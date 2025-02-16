!
!-----------------------------------------------------------------------
!
!  Timing variables
!
!-----------------------------------------------------------------------
!
  REAL ::                                                               &
       cpu0,                                                            &
       cpu_advect, cpu_mixing,   cpu_smlstep, cpu_main,                 &
       cpu_cloud,  cpu_boundary, cpu_copy,    cpu_filter,               &
       cpu_usrout, cpu_init0,    cpu_regrid,  cpu_bndexch,              &
       cpu_bndcint,cpu_update,   cpu_pack,    cpu_unpack,               &
       cpu_simulation
  COMMON /cpu_stat/                                                     &
       cpu_advect, cpu_mixing,   cpu_smlstep, cpu_main,                 &
       cpu_cloud,  cpu_boundary, cpu_copy,    cpu_filter,               &
       cpu_usrout, cpu_init0,    cpu_regrid,  cpu_bndexch,              &
       cpu_bndcint,cpu_update,   cpu_pack,    cpu_unpack,               &
       cpu_simulation
