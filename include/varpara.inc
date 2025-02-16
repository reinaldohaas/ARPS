!
!-----------------------------------------------------------------------
!
!  3DVAR control parameters
!  Input via namelist in INITADAS
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: maxpass = 8
  INTEGER            :: maxin(maxpass)
  COMMON /var_cont/maxin

  INTEGER :: ipass_filt(maxpass), vradius_opt
  REAL    :: hradius(maxpass),vradius(maxpass)
  COMMON /var_refil/ipass_filt, hradius, vradius_opt, vradius

  INTEGER :: chk_opt,assim_opt,cntl_var
  COMMON /var_exprt/chk_opt, assim_opt, cntl_var
  INTEGER :: cntl_var_rh


  INTEGER :: div_opt
  REAL    :: wgt_div_h(maxpass),wgt_div_v(maxpass)
  COMMON /var_diverge/div_opt, wgt_div_h, wgt_div_v

  INTEGER :: smth_flag,thermo_opt
  REAL    :: wgt_smth(maxpass),wgt_thermo(maxpass)
  COMMON /var_smth/smth_flag,wgt_smth,thermo_opt,wgt_thermo
