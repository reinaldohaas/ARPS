!
!  here are the all FORTRAN packing routines from
!  NCAR.  they should work on all machines
!
!  note: we are using the slowest version  of these routines.
!  min and max values are determined and the numbers packed
!  between these values. (we're using packaf, see the NCAR routines).
!
!  installed into the interface: 5 MAY 1992, wcs
!
!  first, our local interface (wrlcm, rdlcm) to the NCAR routines
!

SUBROUTINE wrlcm (arr,apr,nw,npack)
  INCLUDE 'agricpu.inc'
!
!
! ON INPUT:
!  ARR(NW)        AN ARRAY OF REAL NUMBERS
!  NW             ITS LENGTH
!
! ON OUTPUT:
!  APR(1)         MINIMUM OF ARR.
!  APR(2)         THE SCALE FACTOR=S(NPACK-1)/(2*MAX(ABS(ARR)))
!  APR(.)         THE PACKED DATA
!
  DIMENSION arr(nw),apr(*)

!  print*,' npack in wrlcm =' , npack

  IF( npack < 1 .OR. npack > 4 ) THEN
    PRINT*,' wrong value of npack'
    STOP
  END IF
!
  IF (npack == 1) THEN
    DO i=1,nw
      apr(i) = arr(i)
    END DO
  ELSE
    cpu0 = f_cputime()
    CALL packaf( arr,nw,apr,npack,sca1,sca2 )
    cpu_pack = cpu_pack + f_cputime() - cpu0
  END IF
  RETURN
END SUBROUTINE wrlcm
!

SUBROUTINE rdlcm (arr,apr,nw,npack)
  INCLUDE 'agricpu.inc'
!
!
! ON INPUT:
!  APR(.)         THE PACKED ARRAY PRODUCED BY WRLCM
!  NW             THE NUMBER OF ELEMENTS WANTED
!
! ON EXIT:
!  ARR(NW)        THE ARRAY OF REAL NUMBERS
!
!
  DIMENSION arr(nw),apr(*)
!
  IF (npack == 1) THEN
    DO i=1,nw
      arr(i) = apr(i)
    END DO
  ELSE
    cpu0 = f_cputime()
    CALL unpkaf( arr,nw,apr,npack,sca1,sca2 )
    cpu_unpack = cpu_unpack + f_cputime() - cpu0
  END IF
!
  RETURN
END SUBROUTINE rdlcm
!
!  here are the NCAR routines.
!
!***********************************************************************
!
! latest revision       September, 1991.  See "purpose" section
!                    immediately below, and "revision history"
!                    comments at the end of these routines for
!                    more details.
!
! purpose               This version of the packing/unpacking routines
!                    was written in Fortran 77, to be used on CRAY
!                    machines, including XMP, YMP, and CRAY-2 series.
!                    It is intended that this version will replace
!                    the mixed Fortran/Cal version that has been
!                    traditional on NCAR's XMP/YMP machines.
!
!                    This all-Fortran version produces slightly
!                    different results than the FORTRAN/CAL version
!                    unless rounding is turned off on the YMP compiler
!                    cf77 (option"r").  This option is not available
!                    on CRAY-2 machines, which is the reason that the
!                    two sets of routines produce different results.
!
!                    Historically, the purpose of these routines is
!                    to reduce the data volume and possibly the core
!                    requirements associated with a user program by
!                    packing/unpacking several values per Cray word.
!                    See the "algorithm" section below for details.
!
! usage                 call packaf (arn,nrn,apn,npk,amn,amx)
!                    call packbf (arn,nrn,apn,npk,amn,amx)
!                    call packcf (arn,nrn,apn,npk,amn,amx)
!                    call unpkaf (arn,nrn,apn,npk,amn,amx)
!                    call unpkbf (arn,nrn,apn,npk,amn,amx)
!                    call unpkcf (arn,nrn,apn,npk,amn,amx)
!
! arguments             All of the packing and unpacking routines
!                    have the same calling sequence, with
!                    arguments arn, nrn, apn, npk, amn, and amx,
!                    the last two of which may be omitted
!                    except in calls to packcf.  The arguments
!                    may be described as follows:
!
! on input              arn(i),i=1,nrn
!
!                    an array of real numbers - input to the
!                    packers, output from the unpacker
!
!                    nrn     the length of arn - input
!
!                    apn(1)    minimum value - amn
!                    apn(2)    scale factor -
!                              (2**(64/npk)-1)/(amx-amn)
!
!
!                    npk
!
!                    packing density - number of data to be packed
!                    in one output word (2, 3, or 4) - input
!
!                    amn
!
!                    minimum - input to packcf, output from others
!                    (may be omitted from all except packcf calls)
!
!                    amx
!
!                    maximum - input to packcf, output from others
!                    (may be omitted from all except packcf calls)
!
! on output             apn(j) , j=3,4,...n
!
!                    the packed data - npk integers per word -
!                    n = 2+(nrn+npk-1)/npk - output from the
!                    packers, input to the unpacker
!
! special conditions    There are three packing routines to
!                    choose from - packaf, packbf, and packcf.
!                    They differ only in the way in which
!                    amn and amx are computed.
!
!                    Packaf is the safest, but the slowest.
!                    amn and amx are computed by the routine
!                    itself.  Initially, amn is the actual
!                    minimum and amx is the actual maximum.
!                    then, if amn = amx, the values are
!                    adjusted outwards to avoid division by zero.
!                    further, if amn amx differ in sign,
!                    amx is adjusted outwards a little to
!                    ensure that zeroes in the original array
!                    will be returned by the unpacker as zeroes
!                    (mathematically, at least - hardware
!                    round-off makes it impossible to ensure
!                    the return of exact zeroes without
!                    seriously degrading the speed of the
!                    routine - however, the value returned should
!                    never be greater than 1.e-14*abs(amn) ).
!
!                    Packbf is the intermediate choice, faster
!                    but not as safe. amn and amx are computed
!                    by the routine itself.
!                    Initially, amx is set equal to the maximum
!                    absolute value and amn to the negative
!                    of that.  Then, amn and amx are adjusted
!                    as for packaf.  This routine should probably
!                    be used only for data which is known to have
!                    a distribution which is roughly symmetric
!                    around zero.
!
!                    Packcf is the fastest, but most dangerous,
!                    of the three. the user must supply the
!                    initial values of amn and amx.  These values
!                    are not adjusted in any way.
!                    No test is made to ensure that data
!                    values lie between amn and amx - errors may
!                    result.
!
!                    Routines 'unpkaf', 'unpkbf', and 'unpkcf'
!                    are used for the unpacking process.
!                    Actually, they are all synonymns, provided
!                    for the sake of consistency if and when
!                    planned extensions are added to this package.
!
!                    Note that, in any case, these routines
!                    must be used with care.  One wildly
!                    out-of-range datum can destroy the precision
!                    of the rest.  Moreover, data sets with
!                    unlike ranges should be packed separately.
!
! i/o                   none
!
! precision             see algorithm discussion below
!
! required library      none
! files
!
! language              fortran
!
! history               1979: Written by Vince Wayland of Cray Research,
!                    Inc, and made available to NCAR in the same year.
!                    1991: (September) Re-written entirely in fortran
!                    for use on CRAY-2 machines, and also XMP/YMP
!                    machines if desired.  The two versions produce
!                    different results on the XMP/YMP series, as
!                    explained above in the "purpose" section.
!
! algorithm             Given an array of real numbers
!                    (arn(i), i=1,nrn), a packing density
!                    npk (with value 2, 3, or 4), and real numbers
!                    amn and amx (the former less than or equal
!                    to any datum in arn and the latter greater
!                    than or equal to any datum in arn), one can
!                    reduce each element of arn to a positive
!                    integer less than or equal to
!                    2**(64/npk)-1 by means of the
!                    following formula:
!
!                    int(i) = ifix( (arn(i)-amn) *
!                             (2**(64/npk)-1) / (amx-amn) + .5 )
!
!                    the numbers amn and (2**(64/npk)-1) / (amx-amn)
!                    may then be stored as apn(1) and apn(2), followed
!                    by the integers (int(i),i=1,nrn), packed npk per
!                    word in apn(3) through apn(2+(nrn+npk-1)/npk)).
!                    The original (arn(i),i=1,nrn) may subsequently
!                    be reconstructed from the contents of the packed
!                    array apn by means of the following formula:
!
!                    arn(i) = amn + int(i) * (amx-amn) /
!                             (2**(64/npk)-1)
!
!                    Note that the relative error of the reconstructed
!                    values is a function of the packing density.
!                    The error being 1. / (2**(64/npk)-1), i.e.
!                    2.3 e-9, 4.8 e-7, 1.5 e-5 for npk 2, 3, 4
!                    respectively.
!
! portability           not portable.
!
! timing                as of September, 1992 the following timings
!                    obtained on a Cray-YMP8/864.
!
!                    each time given is the time required
!                    to pack or unpack 100,000 numbers, in milliseconds.
!
!                    routine     npk=2      npk=3      npk=4
!                    -------     -----      -----      -----
!
!                    packaf       4.97       5.03       5.11
!                    packbf       3.60       3.67       3.76
!                    packcf       2.06       2.12       2.21
!                    unpkbf       1.79       1.46       1.30
!
!***********************************************************************

SUBROUTINE packaf (arn,nrn,apn,npk,amn,amx)
  DIMENSION arn(nrn),apn(*),po2(3)
  DATA po2 / 4294967295.,2097151.,65535. /
!
! check for certain errors in the arguments.
!
  noa=numarg()
  IF (noa < 4 .OR. noa > 6) GO TO 901
!
  IF (nrn <= 0) GO TO 902
!
  IF (npk < 2 .OR. npk > 4) GO TO 903
!
  xmn=0577777777777777777777B
  xmx=1577777777777777777777B
!
  DO i = 1, nrn
    xmn=MIN(xmn,arn(i))
    xmx=MAX(xmx,arn(i))
  END DO
!
! adjust the values of xmn and xmx - make sure they are not equal.
!
  102 IF (xmn /= xmx) GO TO 104
  xmn=xmn-.5*ABS(xmn)
  xmx=xmx+.5*ABS(xmx)
  IF (xmn /= xmx) GO TO 104
  xmn=-1.
  xmx=+1.
!
! if xmn and xmx differ in sign, further adjust them so that zero will
! be a possible output from the unpacker.
!
  104 IF (xmn*xmx >= 0.) GO TO 105
  IF (ABS(xmn) < ABS(xmx)) THEN
    xmn=xmx-xmx*po2(npk-1)/FLOAT(IFIX( xmx*po2(npk-1)/(xmx-xmn)))
  ELSE
    xmx=xmn-xmn*po2(npk-1)/FLOAT(IFIX(-xmn*po2(npk-1)/(xmx-xmn)))
  END IF
!
! return the computed values of xmn and xmx to the user, if appropriate.
!
  105 IF (noa >= 5) amn=xmn
  IF (noa >= 6) amx=xmx
!
! store the minimum and the scale factor in the output array.
!
  apn(1)=xmn
  apn(2)=po2(npk-1)/(xmx-xmn)
!
! call an appropriate routine to do the actual packing
!
  IF (npk == 2) THEN
    CALL pack2f (arn, nrn, apn)
  ELSE IF (npk == 3) THEN
    CALL pack3f (arn, nrn, apn)
  ELSE
    CALL pack4f (arn, nrn, apn)
  END IF
!
  RETURN
!
! error exits.
!
  901 PRINT 9001
  STOP
!
  902 PRINT 9002
  STOP
!
  903 PRINT 9003
  STOP
!
! formats.
!
  9001 FORMAT (' packaf - wrong number of arguments - stop ')
  9002 FORMAT (' packaf - array length .le. 0 - stop ')
  9003 FORMAT (' packaf - packing density .ne. 2, 3, or 4 - stop ')
!
END SUBROUTINE packaf

SUBROUTINE packbf (arn,nrn,apn,npk,amn,amx)
  DIMENSION arn(nrn),apn(*),po2(3)
  DATA po2 / 4294967295.,2097151.,65535. /
!
! check for certain errors in the arguments.
!
  noa=numarg()
  IF (noa < 4 .OR. noa > 6) GO TO 901
!
  IF (nrn <= 0) GO TO 902
!
!  print*,' npk =', npk

  IF (npk < 2 .OR. npk > 4) GO TO 903
!
  xmx=1577777777777777777777B
!
  DO i = 1, nrn
    xmx=MAX(xmx,arn(i))
  END DO
!
  102 xmn=-xmx
!
! adjust the values of xmn and xmx - make sure they are not equal.
!
  IF (xmn /= xmx) GO TO 103
  xmn=-1.
  xmx=+1.
!
! xmn and xmx differ in sign.  further adjust them so that zero will
! be a possible output from the unpacker.
!
  103 xmx=xmn-xmn*po2(npk-1)/FLOAT(IFIX(-xmn*po2(npk-1)/(xmx-xmn)))
!
! return the computed values of xmn and xmx to the user, if appropriate.
!
  104 IF (noa >= 5) amn=xmn
  IF (noa >= 6) amx=xmx
!
! store the minimum and the scale factor in the output array.
!
  apn(1)=xmn
  apn(2)=po2(npk-1)/(xmx-xmn)
!
! call an appropriate routine to do the actual packing
!
  IF (npk == 2) THEN
    CALL pack2f (arn, nrn, apn)
  ELSE IF (npk == 3) THEN
    CALL pack3f (arn, nrn, apn)
  ELSE
    CALL pack4f (arn, nrn, apn)
  END IF
!
  RETURN
! error exits.
!
  901 PRINT 9001
  STOP
!
  902 PRINT 9002
  STOP
!
  903 PRINT 9003
  STOP
!
! formats.
!
  9001 FORMAT (' packbf - wrong number of arguments - stop ')
  9002 FORMAT (' packbf - array length .le. 0 - stop ')
  9003 FORMAT (' packbf - packing density .ne. 2, 3, or 4 - stop ')
!
END SUBROUTINE packbf

SUBROUTINE packcf (arn,nrn,apn,npk,amn,amx)
!
! this is the subroutine packcf.
!
  DIMENSION arn(nrn),apn(*),po2(3)
  DATA po2 / 4294967295.,2097151.,65535. /
!
! check for certain errors in the arguments.
!
  noa=numarg()
  IF (noa /= 6) GO TO 901
!
  IF (nrn <= 0) GO TO 902
!
  IF (npk < 2 .OR. npk > 4) GO TO 903
!
  IF (amn >= amx) GO TO 904
!
! store the minimum and the scale factor in the output array.
!
  apn(1)=amn
  apn(2)=po2(npk-1)/(amx-amn)
!
! call an appropriate routine to do the actual packing
!
  IF (npk == 2) THEN
    CALL pack2f (arn, nrn, apn)
  ELSE IF (npk == 3) THEN
    CALL pack3f (arn, nrn, apn)
  ELSE
    CALL pack4f (arn, nrn, apn)
  END IF
!
  RETURN
!
! error exits.
!
  901 PRINT 9001
  STOP
!
  902 PRINT 9002
  STOP
!
  903 PRINT 9003
  STOP
!
  904 PRINT 9004
  STOP
!
! formats.
!
  9001 FORMAT (' packcf - wrong number of arguments - stop ')
  9002 FORMAT (' packcf - array length .le. 0 - stop ')
  9003 FORMAT (' packcf - packing density .ne. 2, 3, or 4 - stop ')
  9004 FORMAT (' packcf - minimum .ge. maximum - stop ')
!
END SUBROUTINE packcf

SUBROUTINE pack2f (arn, nrn, apn)
!
! This routine packs 2 floats to the word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 2
  xsf = apn(1)
  xsc = apn(2)
  DO j = 1, npc
    i = j * 2 - 1
    ih = (arn(i  ) - xsf) * xsc + 0.5
    il = (arn(i+1) - xsf) * xsc + 0.5
    apn(j+2) = OR (shiftl (ih, 32), il)
  END DO
!
! Pack up any left-overs
!
  IF (npc * 2 /= nrn) THEN
    j = npc + 1
    i = j * 2 - 1
    ih = (arn(i  ) - xsf) * xsc + 0.5
    apn(j+2) = OR (shiftl (ih, 32), 0)
  END IF
END SUBROUTINE pack2f

SUBROUTINE pack3f (arn, nrn, apn)
!
! This routine packs 3 floats to the word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 3
  xsf = apn(1)
  xsc = apn(2)
  DO j = 1, npc
    i = j * 3 - 2
    ih = (arn(i  ) - xsf) * xsc + 0.5
    im = (arn(i+1) - xsf) * xsc + 0.5
    il = (arn(i+2) - xsf) * xsc + 0.5
    apn(j+2) = OR (shiftl (ih, 42),                                     &
                   OR (shiftl (im, 21), il))
  END DO
!
! Pack up any left-overs
!
  m = nrn - (npc * 3)
  IF (m /= 0) THEN
    j = npc + 1
    i = j * 3 - 2
    ih = (arn(i  ) - xsf) * xsc + 0.5
    IF (m == 1) THEN
      im = 0
    ELSE
      im = (arn(i+1) - xsf) * xsc + 0.5
    END IF
    apn(j+2) = OR (shiftl (ih, 42),                                     &
                   shiftl (im, 21))
  END IF
END SUBROUTINE pack3f

SUBROUTINE pack4f (arn, nrn, apn)
!
! This routine packs 4 floats to the word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 4
  xsf = apn(1)
  xsc = apn(2)
  DO j = 1, npc
    i = j * 4 - 3
    ihh = (arn(i  ) - xsf) * xsc + 0.5
    ihl = (arn(i+1) - xsf) * xsc + 0.5
    ilh = (arn(i+2) - xsf) * xsc + 0.5
    ill = (arn(i+3) - xsf) * xsc + 0.5
    apn(j+2) = OR (shiftl (ihh, 48),                                    &
               OR (shiftl (ihl, 32),                                    &
               OR (shiftl (ilh, 16), ill)))
  END DO
!
! Pack up any left-overs
!
  m = nrn - (npc * 4)
  IF (m > 0) THEN
    j = npc + 1
    i = j * 4 - 3
    ihh = (arn(i  ) - xsf) * xsc + 0.5
    ilh = 0
    ihl = 0
    GO TO (112, 111, 110) M
    110    ilh = (arn(i+2) - xsf) * xsc + 0.5
    111    ihl = (arn(i+1) - xsf) * xsc + 0.5
    112    CONTINUE
    apn(j+2) = OR (shiftl (ihh, 48),                                    &
               OR (shiftl (ihl, 32), shiftl (ilh, 16)))
  END IF
END SUBROUTINE pack4f

SUBROUTINE unpack2f (arn, nrn, apn)
!
! This routine unpacks 2 floats from each word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 2
  xsf = apn(1)
  xsc = 1.0 / apn(2)
  DO j = 1, npc
    i = j * 2 - 1
    ih =      shiftr (apn(j+2), 32)
    il = AND (        apn(j+2), 4294967295)
    arn(i  ) = ih * xsc + xsf
    arn(i+1) = il * xsc + xsf
  END DO
!
! Unpack any left-overs
!
  m = nrn - (npc * 2)
  IF (m > 0) THEN
    j = npc + 1
    i = j * 2 - 1
    ih = shiftr (apn(j+2), 32)
    arn(i) = ih * xsc + xsf
  END IF
END SUBROUTINE unpack2f

SUBROUTINE unpack3f (arn, nrn, apn)
!
! This routine unpacks 3 floats from each word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 3
  xsf = apn(1)
  xsc = 1.0 / apn(2)
  DO j = 1, npc
    i = j * 3 - 2
    ih =      shiftr (apn(j+2), 42)
    im = AND (shiftr (apn(j+2), 21), 2097151)
    il = AND (        apn(j+2),      2097151)
    arn(i  ) = ih * xsc + xsf
    arn(i+1) = im * xsc + xsf
    arn(i+2) = il * xsc + xsf
  END DO
!
! Unpack any left-overs
!
  m = nrn - (npc * 3)
  IF (m > 0) THEN
    j = npc + 1
    i = j * 3 - 2
    ih = shiftr (apn(j+2), 42)
    arn(i) = ih * xsc + xsf
    IF (m > 1) THEN
      im = AND (shiftr (apn(j+2), 21), 2097151)
      arn(i+1) = im * xsc + xsf
    END IF
  END IF
END SUBROUTINE unpack3f

SUBROUTINE unpack4f (arn, nrn, apn)
!
! This routine unpacks 4 floats from each word
!
  DIMENSION arn (nrn), apn (*)
  npc = nrn / 4
  xsf = apn(1)
  xsc = 1.0 / apn(2)
  DO j = 1, npc
    i = j * 4 - 3
    ihh =      shiftr (apn(j+2), 48)
    ihl = AND (shiftr (apn(j+2), 32), 65535)
    ilh = AND (shiftr (apn(j+2), 16), 65535)
    ill = AND (        apn(j+2),      65535)
    arn(i  ) = ihh * xsc + xsf
    arn(i+1) = ihl * xsc + xsf
    arn(i+2) = ilh * xsc + xsf
    arn(i+3) = ill * xsc + xsf
  END DO
!
! Unpack any left-overs
!
  m = nrn - (npc * 4)
  IF (m > 0) THEN
    j = npc + 1
    i = j * 4 - 3
    ihh =      shiftr (apn(j+2), 48)
    ihl = AND (shiftr (apn(j+2), 32), 65535)
    ilh = AND (shiftr (apn(j+2), 16), 65535)
    GO TO (112, 111, 110) M
    110    arn(i+2) = ilh * xsc + xsf
    111    arn(i+1) = ihl * xsc + xsf
    112    arn(i  ) = ihh * xsc + xsf
  END IF
END SUBROUTINE unpack4f

SUBROUTINE unpkbf (arn,nrn,apn,npk,amn,amx)
!
! this is the subroutine unpkbf.
!
  DIMENSION arn(nrn),apn(*),msk(3)
  DATA msk / 4294967295 , 2097151 , 65535 /
!
! 'unpkaf', and 'unpkcf' are all synonyms for 'unpkbf'.
!
  ENTRY unpkaf (arn,nrn,apn,npk,amn,amx)
  ENTRY unpkcf (arn,nrn,apn,npk,amn,amx)
!
! check for errors in the input.
!
  noa=numarg()
  IF (noa < 4 .OR. noa > 6) GO TO 901
!
  IF (nrn <= 0) GO TO 902
!
  IF (npk < 2 .OR. npk > 4) GO TO 903
!
! return the values of amn and amx, if appropriate.
!
  IF (noa >= 5) amn=apn(1)
  IF (noa >= 6) amx=apn(1)+FLOAT(msk(npk-1))/apn(2)
!
! call an appropriate routine to do the actual unpacking
!
  IF (npk == 2) THEN
    CALL unpack2f (arn, nrn, apn)
  ELSE IF (npk == 3) THEN
    CALL unpack3f (arn, nrn, apn)
  ELSE
    CALL unpack4f (arn, nrn, apn)
  END IF
!
  106 RETURN
!
! error exits.
!
  901 PRINT 9001
  STOP
!
  902 PRINT 9002
  STOP
!
  903 PRINT 9003
  STOP
!
! formats.
!
  9001 FORMAT (' unpkf(abc) - wrong number of arguments - stop ')
  9002 FORMAT (' unpkf(abc) - array length .le. 0 - stop ')
  9003 FORMAT (' unpkf(abc) - packing density .ne. 2, 3, or 4 - stop ')
!
! revision history -----------------------------------------------------
!
! september, 1991  Rewrote these routines in Fortran, for use on CRAY-2
!               machines.  This version does not produce results which
!               are identical to results produced by the mixed Fortran
!               /Cal routines used on the YMP machines. The reason is
!               stated above in the "purpose" section of the package
!               documentation.
! ----------------------------------------------------------------------
END SUBROUTINE unpkbf



