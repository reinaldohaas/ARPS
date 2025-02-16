!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                       MP.INC                         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'mp.inc' for ARPS
!
!  This file contains message passing parameters.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Gene Bassett
!  2000/04/27
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Message passing variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: mp_opt          ! Message passing option
                             ! = 0, no message passing
                             ! = 1, use message passing option.

  INTEGER :: nprocs          ! Number of processors.
  INTEGER :: nproc_x         ! Number of processors in x-direction.
  INTEGER :: nproc_y         ! Number of processors in y-direction.
  INTEGER :: nproc_x_in      ! The nproc_x specified in the input file.
  INTEGER :: nproc_y_in      ! The nproc_y specified in the input file.
  INTEGER :: nproc_x_out
  INTEGER :: nproc_y_out
  INTEGER :: loc_x           ! Processor x-location (1 to nproc_x).
  INTEGER :: loc_y           ! Processor y-location (1 to nproc_y).

  INTEGER :: myproc          ! Processor number (0 to nprocs-1).

  INTEGER, PARAMETER :: max_proc=20000  ! Maximum number of processors.

  INTEGER :: proc(max_proc)  ! Processor numbers.

  INTEGER :: max_fopen       ! Maximum number of files allowed open.
  INTEGER :: gentag

  INTEGER, PARAMETER :: FINDX_H = 1, FINDX_R = 2, FINDX_B = 3,          &
                        FINDX_S = 4, FINDX_T = 5, FINDX_A = 6,          &
                        FINDX_I = 7, FINDX_NUM = 7
                        ! FINDX_NUM cannot exceed 10.

  INTEGER :: joindmp(FINDX_NUM)     ! History dump format
                                    ! = 0, dump file for each processor
                                    ! = 1, dump one joined file.

  INTEGER :: readsplit(FINDX_NUM)   ! External data file read option
                                    ! = 0, each processor do its own read
                                    ! = 1, do split on-the-fly

  INTEGER :: readstride      ! = nprocs     if readsplit == 1
                             ! = max_fopen  otherwise
  INTEGER :: dumpstride      ! = nprocs     if joindmp == 1
                             ! = max_fopen  otherwise
  INTEGER :: splitdmp, splitexbc, splitsoil
                             ! Output for smaller pathes than the processor
                             ! subdomain. The size is determined by
                             ! nproc_x_out/nproc_y_out
                             ! splitdmp  - for ARPS history files
                             ! splitexbc - for laternal boundary files
                             ! splitsoil - for soil var files

  COMMON /arpsc005/ mp_opt,nprocs,nproc_x,nproc_y,loc_x,loc_y,          &
                    myproc,proc,max_fopen,gentag,joindmp,readsplit,     &
                    nproc_x_in,nproc_y_in, readstride,     &
                    dumpstride,nproc_x_out,nproc_y_out,                 &
                    splitdmp,splitexbc,splitsoil

  INTEGER, PARAMETER :: tag_w =11,tag_e =12,tag_n =13,tag_s =14
  INTEGER, PARAMETER :: tag_sw=15,tag_se=16,tag_nw=17,tag_ne=18
