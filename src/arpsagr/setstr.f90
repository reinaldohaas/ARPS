SUBROUTINE setstr( mptr )
!
!  this subroutine sets the storage for a grid
!  by allocating space and setting the pointers for the
!  1d, 2d and 3d arrays in common block pntstr
!
  INCLUDE 'agrigrid.inc'
  INCLUDE 'nodal.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
  DIMENSION isp(9)
!
!  get space for all the variables
!
  nx = node(5, mptr)
  ny = node(6, mptr)
  nz = node(14,mptr)
!
!  set the pointers into storage for this grid,
!  first, we'll set the pointers as if they started
!  from storage location zero, after all are set, we
!  save the size (for data dumps, etc...) and allocate the
!  actual storage (if we're not just checking memory size)
!
  IF( verbose4 ) WRITE(6,'('' IN SETSTR, NX,NY,NZ '',3I6)') nx,ny,nz
!
  inptr = 0
  ipint(mptr)    = inptr
  ipreal(mptr)   = ipint(mptr)    + nsint
  ips1d(mptr)    = ipreal(mptr)   + nsreal
  ipx(mptr)      = ips1d(mptr)    + ns1d
  ipy(mptr)      = ipx(mptr)      + nx1d*nx
  ipz(mptr)      = ipy(mptr)      + ny1d*ny
  ipxy(1,mptr)   = ipz(mptr)      + nz1d*nz
  IF( verbose4 ) WRITE(6,'(''  SETSTR POINT 1, IPXY(1,MPTR) = '',I6)')  &
      ipxy(1,mptr)
!
!  we must do a bit more work to figure out the size
!  of the storage for the possibly packed data
!
  IF (nxy2d >= 1) THEN
    DO ivar=2,nxy2d+1
      npk = ipkxy(ivar-1)
      IF ( npk == 1 ) THEN
        ipxy(ivar,mptr) = ipxy(ivar-1,mptr) +                           &
                            nx*ny
      ELSE
        ipxy(ivar,mptr) = ipxy(ivar-1,mptr) +                           &
                               3 + (nx*ny)/npk
      END IF
    END DO
  END IF
  ipxz(1,mptr) = ipxy(nxy2d+1,mptr)
  IF( verbose4 ) WRITE(6,'(''  SETSTR POINT 2, IPXZ(1,MPTR) = '',I6)')  &
      ipxz(1,mptr)
!
  IF (nxz2d >= 1) THEN
    DO ivar=2,nxz2d+1
      npk = ipkxz(ivar-1)
      IF ( npk == 1 ) THEN
        ipxz(ivar,mptr) = ipxz(ivar-1,mptr) +                           &
                            nx*nz
      ELSE
        ipxz(ivar,mptr) = ipxz(ivar-1,mptr) +                           &
                               3 + (nx*nz)/npk
      END IF
    END DO
  END IF
  ipyz(1,mptr) = ipxz(nxz2d+1,mptr)
  IF( verbose4 ) WRITE(6,'(''  SETSTR POINT 3, IPYZ(1,MPTR) = '',I6)')  &
      ipyz(1,mptr)
!
  IF (nyz2d >= 1) THEN
    DO ivar=2,nyz2d+1
      npk = ipkyz(ivar-1)
      IF ( npk == 1 ) THEN
        ipyz(ivar,mptr) = ipyz(ivar-1,mptr) +                           &
                            ny*nz
      ELSE
        ipyz(ivar,mptr) = ipyz(ivar-1,mptr) +                           &
                            3 + (ny*nz)/npk
      END IF
    END DO
  END IF
  ipxyz(1,mptr) = ipyz(nyz2d+1,mptr)
  IF( verbose4 ) WRITE(6,'(''  SETSTR POINT 4, IPXYZ(1,MPTR) = '',I6)') &
      ipxyz(1,mptr)

  IF (nxyz3d >= 1) THEN
    DO ivar=2,nxyz3d+1
      npk = ipkxyz(ivar-1)
      IF ( npk == 1 ) THEN
        ipxyz(ivar,mptr) = ipxyz(ivar-1,mptr) +                         &
                            nx*ny*nz
      ELSE
        ipxyz(ivar,mptr) = ipxyz(ivar-1,mptr) +                         &
                            (3 + (nx*ny)/npk)*nz
      END IF
      IF( verbose4 ) WRITE(6,'(''  SETSTR, POINT 5.1, IVAR, IP= '',2I9)') &
          ivar,ipxyz(ivar,mptr)
    END DO
  END IF
!
! we now know how much space we need
!
  IF ( mptr /= 1 ) THEN
    needsp = ipxyz(nxyz3d+1,mptr)
  ELSE IF ( lexbc == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Set additional storage for arrays which contain external forced
!  boundary conditions for the base grid.
!
!-----------------------------------------------------------------------
!
    ipexbc(1,mptr) = ipxyz(nxyz3d+1,mptr)
    IF( verbose4 ) WRITE(6,'(''  SETSTR POINT 4, IPEXBC(1,MPTR) = '',I6)') &
        ipexbc(1,mptr)

    IF (nexbc3d >= 1) THEN
      DO ivar=2,nexbc3d+1
        npk = ipkexbc(ivar-1)
        IF ( npk == 1 ) THEN
          ipexbc(ivar,mptr) = ipexbc(ivar-1,mptr) +                     &
                              nx*ny*nz
        ELSE
          ipexbc(ivar,mptr) = ipexbc(ivar-1,mptr) +                     &
                              (3 + (nx*ny)/npk)*nz
        END IF
        IF( verbose4 ) WRITE(6,'(''  SETSTR, POINT 5.1, IVAR, IP= '',2I9)') &
            ivar,ipexbc(ivar,mptr)
      END DO
    END IF
    needsp = ipexbc(nexbc3d+1,mptr)
  ELSE
    needsp = ipxyz(nxyz3d+1,mptr)
  END IF

  inptr = igetsp( needsp )
!
!  reset pointers for starting location in the solution array
!
  ipint(mptr)  = ipint(mptr)  + inptr
  ipreal(mptr) = ipreal(mptr) + inptr
  ips1d(mptr)  = ips1d(mptr)  + inptr
  ipx(mptr)    = ipx(mptr)    + inptr
  ipy(mptr)    = ipy(mptr)    + inptr
  ipz(mptr)    = ipz(mptr)    + inptr

  DO i=1,nxy2d+1
    ipxy(i,mptr) = ipxy(i,mptr) + inptr
    IF((ipkxy(i) > 1) .AND. (i <= nxy2d)) THEN
      ispace = 3 + (nx*ny)/ipkxy(i)
      CALL setpck( ipxy(i,mptr),ispace )
    END IF
  END DO

  DO i=1,nxz2d+1
    ipxz(i,mptr) = ipxz(i,mptr) + inptr
    IF((ipkxz(i) > 1) .AND. (i <= nxz2d)) THEN
      ispace = 3 + (nx*nz)/ipkxz(i)
      CALL setpck( ipxz(i,mptr),ispace )
    END IF
  END DO

  DO i=1,nyz2d+1
    ipyz(i,mptr) = ipyz(i,mptr) + inptr
    IF((ipkyz(i) > 1) .AND. (i <= nyz2d)) THEN
      ispace = 3 + (ny*nz)/ipkyz(i)
      CALL setpck( ipyz(i,mptr),ispace )
    END IF
  END DO

  DO i=1,nxyz3d+1
    ipxyz(i,mptr) = ipxyz(i,mptr) + inptr
    IF((ipkxyz(i) > 1) .AND. (i <= nxyz3d)) THEN
      ispace= 3 + (nx*ny)/ipkxyz(i)
      DO k=1,nz
        ink = ipxyz(i,mptr) + (k-1)*ispace
        CALL setpck( ink,ispace )
      END DO
    END IF
  END DO

  IF ( mptr == 1 .AND. lexbc == 1 ) THEN
    DO i=1,nexbc3d+1
      ipexbc(i,mptr) = ipexbc(i,mptr) + inptr
      IF ( ipkexbc(i) > 1 .AND. i <= nexbc3d ) THEN
        ispace = 3 + (nx*ny)/ipkexbc(i)
        DO k=1,nz
          ink = ipexbc(k,mptr) + (k-1)*ispace
          CALL setpck( ink,ispace )
        END DO
      END IF
    END DO
  END IF

!
!  we're finished here
!
  RETURN
END SUBROUTINE setstr
!
!

SUBROUTINE setpck( in,NUMBER )
  INCLUDE 'agrialloc.inc'
!
!  here we set the data in the array a for packing,
!  this allows us to use the igt* functions with packed data
!  before we actually put anything into the packed storage.
!  if we don't do tis the packing routines often croak
!
  a(in) = 100.
  a(in+1) = 100.
  DO i=2,NUMBER-1
    a(in+i) = 0.
  END DO
  RETURN
END SUBROUTINE setpck
