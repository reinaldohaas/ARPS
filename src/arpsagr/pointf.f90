!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtint( mptr, ivar )
!
  INTEGER :: mptr,ivar
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the integer constants array
!  ivar = point to i'th variable, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF(ivar > nsint) THEN
    WRITE(6,*) ' error in igtint, ivar = ',ivar,                        &
               ' while nsint = ',nsint
    STOP
  END IF
  igtint = ipint(mptr) + MAX0(ivar,1) - 1
  RETURN
  END FUNCTION igtint
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtrel( mptr, ivar )
!
  INTEGER :: mptr,ivar
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the real constants array
!  ivar = point to i'th variable, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF(ivar > nsreal) THEN
    WRITE(6,*) ' error in igtrel, ivar = ',ivar,                        &
               ' while nsreal = ',nsreal
    STOP
  END IF
  igtrel = ipreal(mptr) + MAX0(ivar,1) - 1
  RETURN
  END FUNCTION igtrel
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtns1( mptr, iptr, ilen )
!
  INTEGER :: mptr,ivar,ilen
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the real constants array
!  ivar = point to i'th variable, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF( iptr+ilen > ns1d ) THEN
    WRITE(6,*) ' error in igtns1, iptr+ilen = ',iptr+ilen,              &
               ' while ns1d = ',ns1d
    STOP
  END IF
  igtns1 = ips1d(mptr) + MAX0(iptr,1) - 1
  RETURN
  END FUNCTION igtns1
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtnx1( mptr, ivar )
!
  INTEGER :: mptr,ivar
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the 1-D in x storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF(ivar > nx1d) THEN
    WRITE(6,*) ' error in igtnx1, ivar = ',ivar,                        &
               ' while nx1d = ',nx1d
    STOP
  END IF
  igtnx1 = ipx(mptr) + (MAX0(ivar,1) - 1)*node(5,mptr)
  RETURN
  END FUNCTION igtnx1
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtny1( mptr, ivar )
!
  INTEGER :: mptr,ivar
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the 1-D in y storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF(ivar > ny1d) THEN
    WRITE(6,*) ' error in igtny1, ivar = ',ivar,                        &
               ' while ny1d = ',ny1d
    STOP
  END IF
  igtny1 = ipy(mptr) + (MAX0(ivar,1) - 1)*node(6,mptr)
  RETURN
  END FUNCTION igtny1
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtnz1( mptr, ivar )
!
  INTEGER :: mptr,ivar
  INCLUDE 'nodal.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
!
!  return pointer to the 1-D in z storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  we start by checking whether value ivar is possible
!
  IF(ivar > nz1d) THEN
    WRITE(6,*) ' error in igtnz1, ivar = ',ivar,                        &
               ' while nz1d = ',nz1d
    STOP
  END IF
  igtnz1 = ipz(mptr) + (MAX0(ivar,1) - 1)*node(14,mptr)
  RETURN
  END FUNCTION igtnz1
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtnxy( mptr, ivr, inum )
!
  INTEGER :: mptr,ivar,inum
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: pck,nopck
!
!  return pointer to the 2-D in xy storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  inum = number of arrays wanted.
!  note: we may have to unpack some of this stuff.  The present
!  protocol is to see if any need unpacking.  If none need unpacking
!  we just pass back the pointer to the permanent storage location.
!  if some need unpacking we get temp storage for all of them and then
!  unpack as many as necessary.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in igtnxy, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxy2d) THEN
    WRITE(6,*) ' error in igtnxy, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nxy2d = ',nxy2d
    STOP
  END IF
!
!  now see if any are packed or if all are not packed
!
  pck   = .false.
  nopck = .false.
  nx = node(5,mptr)
  ny = node(6,mptr)
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    IF ( ipkxy(i) == 1 ) THEN
      nopck   = .true.
    ELSE IF( ipkxy(i) > 1 .AND. ipkxy(i) <= 4 ) THEN
      pck = .true.
    ELSE
      PRINT*,'Wrong value of ipk. Job stopped in igtnxy.'
      PRINT*,' ipk=', ipkxy(i), i
      STOP
    END IF
  END DO
!
!  now we do as needed to send pointer back.  First, if
!  nopacking, just send back the appropriate pointer
!
  IF ( nopck .AND. .NOT. pck ) THEN
    igtnxy = ipxy(ivar,mptr)
    RETURN
  END IF
!
!  here we must unpack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and unpacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  needsp = nx*ny*inum
  iptr = igetsp( needsp )
  DO i=ivar,ivar+inum-1
    npk = ipkxy(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retxyz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    iloc = (i-ivar)*nx*ny + iptr
    IF ( npk /= 1 ) THEN
      CALL rdlcm( a(iloc),a(ipxy(i,mptr)),nxy,npk )
    ELSE IF( npk > 1 .AND. npk <= 4 ) THEN
      DO ij=1,nxy
        a(ij+iloc-1) = a(ij+ipxy(i,mptr)-1)
      END DO
    END IF
  END DO
  igtnxy = iptr
!
! we're finished here
!
  RETURN
  END FUNCTION igtnxy
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtnxz( mptr, ivr, inum )
!
  INTEGER :: mptr,ivar,inum
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: pck,nopck
!
!  return pointer to the 2-D in xz storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  inum = number of arrays wanted.
!  note: we may have to unpack some of this stuff.  The present
!  protocol is to see if any need unpacking.  If none need unpacking
!  we just pass back the pointer to the permanent storage location.
!  if some need unpacking we get temp storage for all of them and then
!  unpack as many as necessary.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in igtnxz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxz2d) THEN
    WRITE(6,*) ' error in igtnxz, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nxz2d = ',nxz2d
    STOP
  END IF
!
!  now see if any are packed or if all are not packed
!
  pck   = .false.
  nopck = .false.
  nx = node( 5,mptr)
  nz = node(14,mptr)
  nxz = nx*nz
  DO i=ivar,ivar+inum-1
    IF ( ipkxz(i) == 1 ) THEN
      nopck   = .true.
    ELSE IF( ipkxz(i) > 1 .AND. ipkxz(i) <= 4 ) THEN
      pck = .true.
    ELSE
      PRINT*,'Wrong value of ipk. Job stopped in igtnxz.'
      PRINT*,' ipk=', ipkxz(i), i
      STOP
    END IF
  END DO
!
!  now we do as needed to send pointer back.  First, if
!  nopacking, just send back the appropriate pointer
!
  IF ( nopck .AND. .NOT. pck ) THEN
    igtnxz = ipxz(ivar,mptr)
    RETURN
  END IF
!
!  here we must unpack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and unpacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  needsp = nx*nz*inum
  iptr = igetsp( needsp )
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nx*nz + iptr
    npk = ipkxz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retnxz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL rdlcm( a(iloc),a(ipxz(i,mptr)),nxz,npk )
    ELSE IF( npk > 1 .AND. npk <= 4 ) THEN
      DO ij=1,nxz
        a(ij+iloc-1) = a(ij+ipxz(i,mptr)-1)
      END DO
    END IF
  END DO
  igtnxz = iptr
!
! we're finished here
!
  RETURN
  END FUNCTION igtnxz
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtnyz( mptr, ivr, inum )
!
  INTEGER :: mptr,ivar,inum
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: pck,nopck
!
!  return pointer to the 2-D in yz storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  inum = number of arrays wanted.
!  note: we may have to unpack some of this stuff.  The present
!  protocol is to see if any need unpacking.  If none need unpacking
!  we just pass back the pointer to the permanent storage location.
!  if some need unpacking we get temp storage for all of them and then
!  unpack as many as necessary.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in igtnyz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nyz2d) THEN
    WRITE(6,*) ' error in igtnyz, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nyz2d = ',nyz2d
    STOP
  END IF
!
!  now see if any are packed or if all are not packed
!
  pck   = .false.
  nopck = .false.
  ny = node( 6,mptr)
  nz = node(14,mptr)
  nyz = ny*nz
  DO i=ivar,ivar+inum-1
    IF ( ipkyz(i) == 1 ) THEN
      nopck   = .true.
    ELSE IF( ipkyz(i) > 1 .AND. ipkyz(i) <= 4 ) THEN
      pck = .true.
    ELSE
      PRINT*,'Wrong value of ipk. Job stopped in igtnyz.'
      PRINT*,' ipk=', ipkyz(i), i
      STOP
    END IF
  END DO
!
!  now we do as needed to send pointer back.  First, if
!  nopacking, just send back the appropriate pointer
!
  IF ( nopck .AND. .NOT. pck ) THEN
    igtnyz = ipyz(ivar,mptr)
    RETURN
  END IF
!
!  here we must unpack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and unpacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  needsp = ny*nz*inum
  iptr = igetsp( needsp )
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*ny*nz + iptr
    npk = ipkyz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retnyz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL rdlcm( a(iloc),a(ipyz(i,mptr)),nyz,npk )
    ELSE IF( npk > 1 .AND. npk <= 4 ) THEN
      DO ij=1,nyz
        a(ij+iloc-1) = a(ij+ipyz(i,mptr)-1)
      END DO
    END IF
  END DO
  igtnyz = iptr
!
! we're finished here
!
  RETURN
  END FUNCTION igtnyz
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtxyz( mptr, ivr, inum )
!
  INTEGER :: mptr,ivar,inum
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: pck,nopck
!
!  return pointer to the 3-D in xyz storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  inum = number of arrays wanted.
!  note: we may have to unpack some of this stuff.  The present
!  protocol is to see if any need unpacking.  If none need unpacking
!  we just pass back the pointer to the permanent storage location.
!  if some need unpacking we get temp storage for all of them and then
!  unpack as many as necessary.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in igtxyz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxyz3d) THEN
    WRITE(6,*) ' error in igtxyz, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nxyz3d = ',nxyz3d
    STOP
  END IF
!
!  now see if any are packed or if all are not packed
!
  pck   = .false.
  nopck = .false.
  nx = node( 5,mptr)
  ny = node( 6,mptr)
  nz = node(14,mptr)
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    IF ( ipkxyz(i) == 1 ) THEN
      nopck   = .true.
    ELSE IF( ipkxyz(i) > 1 .AND. ipkxyz(i) <= 4 ) THEN
      pck = .true.
    ELSE
      PRINT*,'Wrong value of ipk. Job stopped in igtxyz.'
      PRINT*,' ipk=', ipkxyz(i), i
      STOP
    END IF
  END DO
!
!  now we do as needed to send pointer back.  First, if
!  nopacking, just send back the appropriate pointer
!
  IF ( nopck .AND. .NOT. pck ) THEN
    igtxyz = ipxyz(ivar,mptr)
    RETURN
  END IF
!
!  here we must unpack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and unpacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors,
!  and here it may be particularily  for there's alot of work
!  and we should be able to spread it out over the processors.
!  all the loops in the "do 30" below are independent.
!
  needsp = nxy*nz*inum
  iptr = igetsp( needsp )
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nxy*nz + iptr
    npk = ipkxyz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retxyz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
!
!fpp$ cncall
!
      DO k=1,nz
        ilocout = iloc+nxy*(k-1)
        ilocin  = ipxyz(i,mptr)                                         &
                  + (k-1)*((nxy/npk)+3)
        CALL rdlcm( a(ilocout),a(ilocin),nxy,npk )
      END DO
    ELSE

      DO ijk=1,nxy*nz
        a(ijk+iloc-1) = a(ijk+ipxyz(i,mptr)-1)
      END DO
    END IF

  END DO
  igtxyz = iptr
!
! we're finished here
!
  RETURN
  END FUNCTION igtxyz
!
! --------------------------------------------------------------------
!

  INTEGER FUNCTION igtexbc( mptr, ivr, inum )

  IMPLICIT NONE

  INTEGER :: mptr,ivr,inum
  INTEGER :: nx,ny,nz,nxy
  INTEGER :: npk,iptr,iloc,ilocin,ilocout
  INTEGER :: i,j,k,ijk,ivar

  INTEGER :: needsp
  INTEGER :: igetsp

  LOGICAL :: pck,nopck

  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
!
!  return pointer to the 3-D in xyz storage arrays
!  ivar = point to i'th array, if zero, point to the first.
!  inum = number of arrays wanted.
!  note: we may have to unpack some of this stuff.  The present
!  protocol is to see if any need unpacking.  If none need unpacking
!  we just pass back the pointer to the permanent storage location.
!  if some need unpacking we get temp storage for all of them and then
!  unpack as many as necessary.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in igtexbc, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nexbc3d) THEN
    WRITE(6,*) ' error in igtexbc, ivar = ',ivar,                       &
               ' inum = ',inum,                                         &
               ' while nexbc3d = ',nexbc3d
    STOP
  END IF
!
!  now see if any are packed or if all are not packed
!
  pck   = .false.
  nopck = .false.
  nx = node( 5,mptr)
  ny = node( 6,mptr)
  nz = node(14,mptr)
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    IF ( ipkexbc(i) == 1 ) THEN
      nopck   = .true.
    ELSE IF( ipkexbc(i) > 1 .AND. ipkexbc(i) <= 4 ) THEN
      pck = .true.
    ELSE
      PRINT*,'Wrong value of ipk. Job stopped in igtexbc.'
      PRINT*,' ipk=', ipkexbc(i), i
      STOP
    END IF
  END DO
!
!  now we do as needed to send pointer back.  First, if
!  nopacking, just send back the appropriate pointer
!
  IF ( nopck .AND. .NOT. pck ) THEN
    igtexbc = ipexbc(ivar,mptr)
    RETURN
  END IF
!
!  here we must unpack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and unpacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors,
!  and here it may be particularily  for there's alot of work
!  and we should be able to spread it out over the processors.
!  all the loops in the "do 30" below are independent.
!
  needsp = nxy*nz*inum
  iptr = igetsp( needsp )
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nxy*nz + iptr
    npk = ipkexbc(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in igtexbc.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
!
!fpp$ cncall
!
      DO k=1,nz
        ilocout = iloc+nxy*(k-1)
        ilocin  = ipexbc(i,mptr)                                        &
                  + (k-1)*((nxy/npk)+3)
        CALL rdlcm( a(ilocout),a(ilocin),nxy,npk )
      END DO
    ELSE

      DO ijk=1,nxy*nz
        a(ijk+iloc-1) = a(ijk+ipexbc(i,mptr)-1)
      END DO
    END IF

  END DO
  igtexbc = iptr
!
! we're finished here
!
  RETURN
  END FUNCTION igtexbc
!
! --------------------------------------------------------------------
!

SUBROUTINE retnxy( mptr, ivr, inum, iptr, resetd )
!
  INTEGER :: mptr,ivar,inum,iptr
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: resetd
!
!  this routine places data back into the permanent storage
!  area if the data was previously unpacked into temp storage
!  and if resetd is true.
!
!  if reset d is false it just recovers the temp space, if any was used
!
!  ivar = point to "ivar" xy 2-d array
!  inum = number of arrays stored at iptr beginning with ivar.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in retnxy, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxy2d) THEN
    WRITE(6,*) ' error in retnxy, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nxy2d = ',nxy2d
    STOP
  END IF
!
!  simple check to see if temp space is used,
!  if not, then just return
!
  IF( iptr == ipxy(ivar,mptr) ) RETURN
!
!  at this point we can return the temp space for it
!  won't be used before we leave
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  isp = nx*ny*inum
  CALL reclam( iptr,isp )
!
!  return if we don't need to reset the permanent storage
!
  IF( .NOT. resetd ) RETURN
!
!  here we must repack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and repacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nx*ny + iptr
    npk = ipkxy(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retnxy.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL fillrw( a(iloc),nx,ny,idmxy(1,i),idmxy(2,i),1 )
      CALL wrlcm( a(iloc),a(ipxy(i,mptr)),nxy,npk )
    ELSE
      DO ij=1,nxy
        a(ij+ipxy(i,mptr)-1) = a(ij+iloc-1)
      END DO
    END IF
  END DO
!
! we're finished here
!
  RETURN
END SUBROUTINE retnxy
!
! --------------------------------------------------------------------
!

SUBROUTINE retnxz( mptr, ivr, inum, iptr, resetd )
!
  INTEGER :: mptr,ivar,inum,iptr
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: resetd
!
!  this routine places data back into the permanent storage
!  area if the data was previously unpacked into temp storage
!  and if resetd is true.
!
!  if reset d is false it just recovers the temp space, if any was used
!
!  ivar = point to "ivar" xz 2-d array
!  inum = number of arrays stored at iptr beginning with ivar.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in retnxz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxz2d) THEN
    WRITE(6,*) ' error in retnxz, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nxz2d = ',nxz2d
    STOP
  END IF
!
!  simple check to see if temp space is used,
!  if not, then just return
!
  IF( iptr == ipxz(ivar,mptr) ) RETURN
!
!  at this point, we can return the temp space, for it
!  won't be used before we leave
!
  nx = node( 5,mptr)
  nz = node(14,mptr)
  isp = nx*nz*inum
  CALL reclam( iptr,isp )
!
!  return if we don't need to reset the permanent storage
!
  IF( .NOT. resetd ) RETURN
!
!  here we must repack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and repacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  nxz = nx*nz
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nx*nz + iptr
    npk = ipkxz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retnxz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL fillrw( a(iloc),nx,nz,idmxz(1,i),idmxz(2,i),1 )
      CALL wrlcm( a(iloc),a(ipxz(i,mptr)),nxz,npk )
    ELSE
      DO ij=1,nxz
        a(ij+ipxz(i,mptr)-1) = a(ij+iloc-1)
      END DO
    END IF
  END DO
!
! we're finished here
!
  RETURN
END SUBROUTINE retnxz
!
! --------------------------------------------------------------------
!

SUBROUTINE retnyz( mptr, ivr, inum, iptr, resetd )
!
  INTEGER :: mptr,ivar,inum,iptr
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: resetd
!
!  this routine places data back into the permanent storage
!  area if the data was previously unpacked into temp storage
!  and if resetd is true.
!
!  if reset d is false it just recovers the temp space, if any was used
!
!  ivar = point to "ivar" yz 2-d array
!  inum = number of arrays stored at iptr beginning with ivar.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in retnyz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nyz2d) THEN
    WRITE(6,*) ' error in retnyz, ivar = ',ivar,                        &
               ' inum = ',inum,                                         &
               ' while nyz2d = ',nyz2d
    STOP
  END IF
!
!  simple check to see if temp space is used,
!  if not, then just return
!
  IF( iptr == ipyz(ivar,mptr) ) RETURN
!
!  at this point, we can return the temp space, for it
!  won't be used before we leave
!
  ny = node( 6,mptr)
  nz = node(14,mptr)
  isp = ny*nz*inum
  CALL reclam( iptr,isp )
!
!  return if we don't need to reset the permanent storage
!
  IF( .NOT. resetd ) RETURN
!
!  here we must repack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and repacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  nyz = ny*nz
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*ny*nz + iptr
    npk = ipkyz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retnyz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL fillrw( a(iloc),ny,nz,idmyz(1,i),idmyz(2,i),1 )
      CALL wrlcm( a(iloc),a(ipyz(i,mptr)),nyz,npk )
    ELSE
      DO ij=1,nyz
        a(ij+ipyz(i,mptr)-1) = a(ij+iloc-1)
      END DO
    END IF
  END DO
!
! we're finished here
!
  RETURN
END SUBROUTINE retnyz
!
! --------------------------------------------------------------------
!

SUBROUTINE retxyz( mptr, ivr, inum, iptr, resetd )
!
  INTEGER :: mptr,ivar,inum,iptr
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: resetd
!
!  this routine places data back into the permanent storage
!  area if the data was previously unpacked into temp storage
!  and if resetd is true.
!
!  if resetd is false it just recovers the temp space, if any was used
!
!  ivar = point to "ivar" xyz 3-d array
!  inum = number of arrays stored at iptr beginning with ivar.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in retxyz, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nxyz3d) THEN
    WRITE(6,*) ' error in retnxyz, ivar = ',ivar,                       &
               ' inum = ',inum,                                         &
               ' while nxyz3d = ',nxyz3d
    STOP
  END IF
!
!  simple check to see if temp space is used,
!  if not, then just return
!
  IF( iptr == ipxyz(ivar,mptr) ) RETURN
!
!  at this point, we can return the temp space, for it
!  won't be used before we leave
!
  nx = node( 5,mptr)
  ny = node( 6,mptr)
  nz = node(14,mptr)
  isp = nx*ny*nz*inum
  IF (.true.) THEN
    PRINT *, 'calling reclam from within retxyz...'
    PRINT *, ' variable number = ', ivr
    PRINT *, 'iptr , space = ', iptr, isp
  END IF

  CALL reclam( iptr,isp )
!
!  return if we don't need to reset the permanent storage
!
  IF( .NOT. resetd ) RETURN
!
!  here we must repack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and repacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nx*ny*nz + iptr
    npk = ipkxyz(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retxyz.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL fillrw( a(iloc),nx,ny,idmxyz(1,i),idmxyz(2,i),nz )
!
!fpp$ cncall
!
      DO k=1,nz
        ilocout = iloc+nx*ny*(k-1)
        ilocin  = ipxyz(i,mptr)                                         &
                  + (k-1)*((nxy/npk)+3)
        CALL wrlcm( a(ilocout),a(ilocin),nxy,npk )
      END DO
    ELSE
      DO ijk=1,nxy*nz
        a(ijk+ipxyz(i,mptr)-1) = a(ijk+iloc-1)
      END DO
    END IF
  END DO
!
! we're finished here
!
  RETURN
END SUBROUTINE retxyz
!
! --------------------------------------------------------------------
!

SUBROUTINE retexbc( mptr, ivr, inum, iptr, resetd )
!
  INTEGER :: mptr,ivar,inum,iptr
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agrialloc.inc'
  LOGICAL :: resetd
!
!  this routine places data back into the permanent storage
!  area if the data was previously unpacked into temp storage
!  and if resetd is true.
!
!  if resetd is false it just recovers the temp space, if any was used
!
!  ivar = point to "ivar" exbc 3-d array
!  inum = number of arrays stored at iptr beginning with ivar.
!
!  we start by checking whether value ivar  and inum are possible
!
  ivar = MAX0(ivr,1)
  IF(inum < 1) THEN
    WRITE(6,*) ' error in retexbc, inum = ',inum
    STOP
  END IF
!
  IF(ivar+inum-1 > nexbc3d) THEN
    WRITE(6,*) ' error in retnexbc, ivar = ',ivar,                      &
               ' inum = ',inum,                                         &
               ' while nexbc3d = ',nexbc3d
    STOP
  END IF
!
!  simple check to see if temp space is used,
!  if not, then just return
!
  IF( iptr == ipexbc(ivar,mptr) ) RETURN
!
!  at this point, we can return the temp space, for it
!  won't be used before we leave
!
  nx = node( 5,mptr)
  ny = node( 6,mptr)
  nz = node(14,mptr)
  isp = nx*ny*nz*inum
  IF (.true.) THEN
    PRINT *, 'calling reclam from within retexbc...'
    PRINT *, ' variable number = ', ivr
    PRINT *, 'iptr , space = ', iptr, isp
  END IF

  CALL reclam( iptr,isp )
!
!  return if we don't need to reset the permanent storage
!
  IF( .NOT. resetd ) RETURN
!
!  here we must repack some of the stuff,
!  this is the slow version that checks on one array at a
!  time and repacks if necessary, else just does a copy.
!  we'll have to see how fast it is on multiple processors
!
  nxy = nx*ny
  DO i=ivar,ivar+inum-1
    iloc = (i-ivar)*nx*ny*nz + iptr
    npk = ipkexbc(i)

    IF( npk < 1 .OR. npk > 4 ) THEN
      PRINT*,'Wrong value of ipk. Job stopped in retexbc.'
      PRINT*,' ipk=', npk , i
      STOP
    END IF

    IF ( npk /= 1 ) THEN
      CALL fillrw( a(iloc),nx,ny,idmexbc(1,i),idmexbc(2,i),nz )
!
!fpp$ cncall
!
      DO k=1,nz
        ilocout = iloc+nx*ny*(k-1)
        ilocin  = ipexbc(i,mptr)                                        &
                  + (k-1)*((nxy/npk)+3)
        CALL wrlcm( a(ilocout),a(ilocin),nxy,npk )
      END DO
    ELSE
      DO ijk=1,nxy*nz
        a(ijk+ipexbc(i,mptr)-1) = a(ijk+iloc-1)
      END DO
    END IF
  END DO
!
! we're finished here
!
  RETURN
END SUBROUTINE retexbc
!
! --------------------------------------------------------------------
!

SUBROUTINE fillrw( a,nx,ny,nxad,nyad,nz )
  INTEGER :: nx,ny,nxad,nyad,nz
  DIMENSION a(nx,ny,nz)
!
!  this routine fills out the outer rows with values before packing.
!  we do this so that the outer rows won't contain any spurious
!  minima or maxima if the outer rows are not used for data.
!  this allows us to retain the maximum precision in the packing and
!  unpacking process.
!
  nxa = nx-nxad
  nya = ny-nyad
  IF( nxa < nx ) THEN
    DO i=nxa+1,nx
      DO k=1,nz
        DO j=1,ny
          a(i,j,k) = a(nxa,j,k)
        END DO
      END DO
    END DO
  END IF
!
  IF( nya < ny ) THEN
    DO j=nya+1,ny
      DO k=1,nz
        DO i=1,nx
          a(i,j,k) = a(i,nya,k)
        END DO
      END DO
    END DO
  END IF
!
  RETURN
END SUBROUTINE fillrw
