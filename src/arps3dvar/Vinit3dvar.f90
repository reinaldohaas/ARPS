!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITADAS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE init3dvar(cntl_var_rh_out,npass, namelist_filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in analysis variables in namelist format from standard input.
!
!  AUTHOR:
!  Jidong GAO, add 3dvar input parameter, 2001
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'varpara.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'
!
  INTEGER :: nx,ny,nz          ! ARPS grid size
  INTEGER :: nstyps            ! Maximum number of soil types per grid point.
  INTEGER :: nt                ! Number of time levels of data
  INTEGER :: i
  INTEGER, INTENT(IN) :: npass
!
!
!-----------------------------------------------------------------------
!
!  3DVAR namelists
!
!-----------------------------------------------------------------------
!
  NAMELIST /var_const/maxin

  NAMELIST /var_refil/ipass_filt,hradius,vradius_opt,vradius

  NAMELIST /var_exprt/chk_opt,assim_opt,cntl_var,cntl_var_rh

  NAMELIST /var_diverge/div_opt,wgt_div_h,wgt_div_v

  NAMELIST /var_smth/smth_flag, wgt_smth

  NAMELIST /var_thermo/thermo_opt, wgt_thermo
!
  INTEGER :: cntl_var_rh_out
  CHARACTER(LEN=*), INTENT(IN) :: namelist_filename

  INTEGER :: unum
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  unum = 0
  IF (myproc == 0) THEN
    IF (LEN_TRIM(namelist_filename) <= 0 .OR. namelist_filename == ' ') THEN
      unum = 5
      WRITE(6,'(2(1x,a,/))') 'Waiting ARPS3DVAR namelist from standard input ... ', &
                             '========================================'
    ELSE
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(namelist_filename),STATUS='OLD',FORM='FORMATTED')
      WRITE(6,'(1x,3a,/,1x,a,/)') 'Reading ARPS3DVAR namelist from file - ',   &
              TRIM(namelist_filename),' ... ','========================================'
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Assign default values to the ADAS input variables
!
!-----------------------------------------------------------------------
!
   maxin(1) = 50
   maxin(2) = 20
   maxin(3) = 30
   ipass_filt(1) = 1
   ipass_filt(2) = 1
   ipass_filt(3) = 1
   hradius(1) = 25.0
   hradius(2) = 25.0
   hradius(3) = 15.0
   vradius_opt= 1
   vradius(1) = 2
   vradius(2) = 2
   vradius(3) = 2
   chk_opt = 0
   assim_opt= 1
   cntl_var = 0
   div_opt = 0
   DO i=1,maxpass
     wgt_div_h(i) = -1.0
     wgt_div_v(i) = -1.0
   ENDDO
   smth_flag  =0
   wgt_smth   = 0.0
   thermo_opt = 0
   wgt_thermo = 0.0
!
!-----------------------------------------------------------------------
!
!  read in ADAS namelists
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ(unum, var_const,  END=350)
    WRITE(6,*) 'Namelist block 3dvar_const sucessfully read.'

    READ(unum, var_refil,  END=350)
    WRITE(6,*) 'Namelist block 3dvar_refil sucessfully read.'

    READ(unum, var_exprt,  END=350)
    WRITE(6,*) 'Namelist block 3dvar_exprt sucessfully read.'

    READ(unum, var_diverge,  END=350)
    WRITE(6,*) 'Namelist block 3dvar_divergence sucessfully read.'

    READ(unum, var_smth,  END=350)
    WRITE(6,*) 'Namelist block 3dvar_smoothness sucessfully read.'
    350 CONTINUE
  END IF
  CALL mpupdatei(maxin,      maxpass)
  CALL mpupdatei(ipass_filt, maxpass)
  CALL mpupdater(hradius,    maxpass)
  CALL mpupdatei(vradius,    maxpass)
  CALL mpupdatei(vradius_opt,  1    )

  CALL mpupdatei(chk_opt,     1)
  CALL mpupdatei(assim_opt,   1)
  CALL mpupdatei(cntl_var,    1)
  CALL mpupdatei(cntl_var_rh, 1)

  CALL mpupdatei(div_opt,   1)
  CALL mpupdater(wgt_div_h, maxpass)
  CALL mpupdater(wgt_div_v, maxpass)

  CALL mpupdatei(smth_flag, 1)
  CALL mpupdater(wgt_smth,  maxpass)

  CALL mpupdatei(thermo_opt, 1)
  CALL mpupdater(wgt_thermo, maxpass)

  DO i=1,npass
    hradius(i)= hradius(i)*1000.0/dx
    IF (myproc == 0) THEN
      WRITE(*,*) 'Analysis pass:', i
      WRITE(*,*) 'The horizontal influence radius:',hradius(i),'grid points'
      WRITE(*,*) 'The vertical influence radius option vradius_opt=',vradius_opt
      IF (vradius_opt == 1 ) THEN
        WRITE(*,*) 'The vertical influence radius:',vradius(i), 'grid points'
      ELSE IF (vradius_opt == 2 ) THEN
        WRITE(*,*) 'The vertical influence radius:',vradius(i), 'km'
      ELSE
        WRITE(*,*) 'vradius_opt = ', vradius_opt, ' is not supported.'
        WRITE(*,*) 'ARPS3DVAR will stop here'
        CALL arpsstop('ERROR: Wrong vradius_opt.',1)
      END IF
    END IF
  END DO

  cntl_var_rh_out=cntl_var_rh

!-----------------------------------------------------------------------
!
  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE(unum)
    CALL retunit(unum)
  END IF

  RETURN
END SUBROUTINE init3dvar
