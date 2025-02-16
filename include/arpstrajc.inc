
  INTEGER :: trajc_plt_opt,ntimes

  INTEGER, PARAMETER :: nmax_times =20
  
  INTEGER, PARAMETER :: npoints_max=10000
  INTEGER, PARAMETER :: ntrajcs_max=30
  INTEGER :: traj_col(ntrajcs_max)
  
  REAL :: xtrajc(npoints_max,ntrajcs_max,nmax_times)
  REAL :: ytrajc(npoints_max,ntrajcs_max,nmax_times)
  REAL :: ztrajc(npoints_max,ntrajcs_max,nmax_times)
  REAL :: ttrajc(npoints_max)
  INTEGER :: ntrajcs(nmax_times), npoints_in(nmax_times)

  INTEGER :: npoints_bgn(ntrajcs_max,nmax_times)
  INTEGER :: npoints_end(ntrajcs_max,nmax_times)

  INTEGER :: itrajc_index,h_follow_trajc

  COMMON /trajc_data/ xtrajc,ytrajc,ztrajc,ttrajc,ntrajcs,npoints_in, &
          npoints_bgn,npoints_end,itrajc_index,h_follow_trajc
  
  INTEGER :: trajc_lbl_opt,trajc_lbl_frq,trajc_mkr_typ,trajc_mkr_frq,trajc_lbl_fmt
  REAL :: trajc_lbl_siz, trajc_mkr_siz
  REAL :: trajc_lbl_number
  
  REAL :: trajc_plt_bgn_time,trajc_plt_end_time
  REAL :: trajc_plotting_startime,trajc_plotting_stoptime
  INTEGER :: ntrajc_start, ntrajc_end, ntrajc_stride
  
  CHARACTER(LEN=256) :: trajc_fn_in(nmax_times)
  COMMON  /trajectopt/ trajc_plt_opt,trajc_plt_bgn_time,trajc_plt_end_time,traj_col, &
           trajc_fn_in,ntimes,trajc_lbl_opt,trajc_lbl_frq, &
           trajc_lbl_siz,trajc_mkr_typ,trajc_mkr_frq,trajc_mkr_siz,trajc_lbl_fmt, &
           ntrajc_start, ntrajc_end, ntrajc_stride

