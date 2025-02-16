!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE module_hdfsubdomain_namelist

  INTEGER, PARAMETER :: MAXHISFILE = 200
  INTEGER, PARAMETER :: MAXDATASET = 100

  INTEGER            :: nhisfile
  CHARACTER(LEN=256) :: hisfile(MAXHISFILE)
  CHARACTER(LEN=256) :: outfile(MAXHISFILE)

  INTEGER            :: nxpnt, nypnt

  CHARACTER(LEN=80)  :: runname

  CHARACTER(LEN=256) :: dirname

  CHARACTER(LEN=20)  :: varnames(MAXDATASET)

  INTEGER            :: nvarout

  INTEGER            :: readyfl

  INTEGER            :: ncompressx, ncompressy
  INTEGER            :: nprocx_out, nprocy_out
  INTEGER            :: nprocx_lw,  nprocy_lw

  LOGICAL            :: inpatch, outpatch

  INTEGER            :: lvldbg

  CONTAINS

  SUBROUTINE read_namelist_params(iamroot,untin, untout, istatus)

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: iamroot
    INTEGER, INTENT(IN)  :: untin, untout
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INCLUDE 'mp.inc'

!-----------------------------------------------------------------------

    INTEGER :: nprocx_in, nprocy_in

    INTEGER :: hdmpinopt, grdbasopt
    CHARACTER(LEN=256) :: hdmpfheader, hdmpftrailer
    REAL               :: tintv_dmpin, tbgn_dmpin, tend_dmpin
    CHARACTER(LEN=256) :: grdbasfn

    NAMELIST /message_passing/ nproc_x, nproc_y


    NAMELIST /history_data/ hdmpinopt, hdmpfheader, hdmpftrailer,       &
                            tintv_dmpin, tbgn_dmpin, tend_dmpin,        &
                            nhisfile, grdbasfn, hisfile,                &
                            nprocx_in,nprocy_in,nprocx_lw, nprocy_lw

    NAMELIST /output_dims/ nxpnt, nypnt

    NAMELIST /jobname/ runname

    NAMELIST /output/ dirname, varnames, grdbasopt, readyfl, lvldbg,  &
                      nprocx_out, nprocy_out

!-----------------------------------------------------------------------

    INTEGER :: nf, nv
    INTEGER :: lstr, lfnkey, n, n1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nprocx_in  = 1
    nprocy_in  = 1

    nprocx_out = 1
    nprocy_out = 1

    nprocx_lw = 1
    nprocy_lw = 1

    nhisfile = 0
    hisfile = ' '

    nxpnt = 0
    nypnt = 0

    runname = 'hdfsubdomain'

    dirname  = './'
    varnames = ' '
    grdbasopt = 1
    nvarout  = 0
    readyfl  = 0
    lvldbg   = 0

    inpatch  = .FALSE.
    outpatch = .FALSE.

!--------------------------- Message Passing ---------------------------
    IF (iamroot) THEN
      READ(untin,NML=message_passing)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block message_passing successfully read.'
      max_fopen = nproc_x*nproc_y
    END IF
    CALL mpupdatei(nproc_x,1)
    CALL mpupdatei(nproc_y,1)

    CALL mpinit_var

!--------------------------- history_data ------------------------------

    IF (iamroot) THEN

      READ(UNTIN,NML=history_data,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block history_data successfully read.'

      IF( hdmpinopt == 1 ) THEN
        CALL gthinfns(hdmpfheader,hdmpftrailer,3,                       &
                  tintv_dmpin, tbgn_dmpin, tend_dmpin,                  &
                  grdbasfn,hisfile,nhisfile)
      END IF

      WRITE(UNTOUT,'(/3x,a,a)')'The grid base-state file name is ',     &
                  TRIM(grdbasfn)

      DO nf=1,nhisfile
        WRITE(UNTOUT,'(3x,a,i5,a,a)')                                   &
        'History file      No. ',nf,'   is ',trim(hisfile(nf))
      END DO

!      nhisfile = nhisfile + 1
!      hisfile(nhisfile) = grdbasfn

!--------------------------- output_dims -------------------------------

      READ(UNTIN,NML=output_dims,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block output_dims successfully read.'

      IF (MOD(nxpnt,2) /= 0 .OR. MOD(nypnt,2) /= 0) THEN
        WRITE(UNTOUT,'(1x,a,/,1x,a,2(I4,a))')                           &
          'ERROR: In order for the output domain has the same center '  &
          //'as the input domain, parameter nxpnt & nypnt must be even' &
          //' numbers.', &
          '       Found nxpnt = ',nxpnt,', nypnt = ',nypnt,'.'
        GO TO 999
      END IF

!--------------------------- jobname ------------------------------

      READ(UNTIN,NML=jobname,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block jobname successfully read.'

 !--------------------------- output ------------------------------

      READ(UNTIN,NML=output,ERR=999)
      WRITE(UNTOUT,'(1x,a)') 'Namelist block output successfully read.'

      lstr = LEN_TRIM(dirname)
      IF (lstr < 1) THEN
        dirname = './'
      ELSE IF (dirname(lstr:lstr) /= '/') THEN
        dirname(lstr+1:lstr+1) = '/'
      END IF

      nvarout = 0
      DO nv = 1, MAXDATASET
        IF (LEN_TRIM(varnames(nv)) > 0) THEN
          nvarout = nvarout + 1
        ELSE
          EXIT
        END IF
      END DO

      IF(grdbasopt /= 0 ) THEN
        nhisfile = nhisfile + 1
        hisfile(nhisfile) = grdbasfn
      END IF

!-----------------------------------------------------------------------
!
! Number of file to be read by each processor
!
!-----------------------------------------------------------------------

      IF (mp_opt == 0) THEN
        readsplit(:) = 0
        ncompressx = nprocx_in
        ncompressy = nprocy_in
      ELSE
        IF (nprocx_in == 1 .AND. nprocy_in == 1) THEN
          readsplit(:) = 1
          ncompressx = 1
          ncompressy = 1
          IF (nprocs > 1) THEN
            WRITE(UNTOUT,'(1x,a,/,8x,a)')                               &
               'ERROR: It has no any benefit for a MPI job with joined file.', &
               'Please use serial version of this program. Program stopped.'
            istatus = -2
          END IF
        ELSE
          IF (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0) THEN
            WRITE(UNTOUT,'(1x,a,/,2(1x,2(a,I4)),/)')                    &
        'ERROR: Number of input patches must be a multiple of the '     &
        //'working number of patches in each dirction.',                &
        '       nprocx_in = ',nprocx_in,', nprocy_in = ',nprocy_in,     &
        '       nproc_x   = ',nproc_x,  ', nproc_y   = ',nproc_y
            GO TO 999
          END IF

          ncompressx = nprocx_in / nproc_x
          ncompressy = nprocy_in / nproc_y
        END IF
      END IF

      IF (nprocx_in  > 1 .OR. nprocy_in  > 1) inpatch  = .TRUE.
      IF (nprocx_out > 1 .OR. nprocy_out > 1) outpatch = .TRUE.

      WRITE(6,*)

      GOTO 900
!-----------------------------------------------------------------------
!
! Error with reading namelist
!
!-----------------------------------------------------------------------
      999 CONTINUE
      istatus = -1

      900 CONTINUE
    END IF
    CALL mpupdatei(istatus,1)

    IF (istatus /= 0) THEN
      IF (IAMROOT) WRITE(UNTOUT,'(1x,a,/,1x,a)') 'Error reading NAMELIST file.', &
         'Job stopped in read_namelist_params.'

      CALL arpsstop('ERROR: Namelist input error',1)
    END IF

!-----------------------------------------------------------------------
!
! Message passing of parameters
!
!-----------------------------------------------------------------------

    CALL mpupdatei(readsplit,FINDX_NUM)

    CALL mpupdatei(nprocx_out,1)
    CALL mpupdatei(nprocy_out,1)

    CALL mpupdatei(ncompressx,1)
    CALL mpupdatei(ncompressy,1)

    CALL mpupdatei(nprocx_lw,1)
    CALL mpupdatei(nprocy_lw,1)

    CALL mpupdatei(nhisfile,1)
    CALL mpupdatec(hisfile,256*nhisfile)

    CALL mpupdatei(nxpnt,1)
    CALL mpupdatei(nypnt,1)

    CALL mpupdatec(runname,80)
    CALL mpupdatec(dirname,256)

    CALL mpupdatei(nvarout,1)
    CALL mpupdatec(varnames,20*nvarout)

    CALL mpupdatei(readyfl,1)

    DO n = 1, nhisfile
      n1 = INDEX(hisfile(n),'.',.TRUE.)
      IF (n1 < 1) n1 = 1
      CALL gtlfnkey( runname, lfnkey )

      WRITE(outfile(n),'(3a)') TRIM(dirname),runname(1:lfnkey),         &
                               TRIM(hisfile(n)(n1:))
    END DO

    RETURN
  END SUBROUTINE read_namelist_params

END MODULE module_hdfsubdomain_namelist
