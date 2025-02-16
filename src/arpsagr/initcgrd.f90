SUBROUTINE domain( nxc,nyc,nzc , timestart)
  INCLUDE 'nodal.inc'
!
!  this routine sets up the coarse grid when the calculation
!  is not being restarted from a previous run.
!
  mstart=nodget(dum)
  node(4,mstart)=0
  lstart(1)=mstart
  lback(1)=mstart

  WRITE(6,'(1x,a)')                                                     &
      'The max. number of nesting levels was set to 10.'

  mptrc = mstart

  rnode(1, mptrc) = 0.
  rnode(2, mptrc) = 0.
  rnode(3, mptrc) = 0.

  node(4, mptrc) = 1
  node(5, mptrc) = nxc  ! base grid nx
  node(6, mptrc) = nyc  ! base grid ny
  node(7, 1)=0
  node(8, 1)=0
  node(10,mstart) = 0
  node(14,mptrc) = nzc  ! base grid nz
!
! not used nodes
!
  node(16,1)=0
  node(17,mptrc) = 2
  node(18,mptrc) = 2
  node(19,mptrc) = 3
  node(20,mptrc) = 3
  node(21,mptrc) = 3
  node(22,mptrc) = 3

  CALL setstr( mstart )
  CALL settmp

  nestgrd0 = 1

  CALL arpsinit( mstart , nestgrd0 )

  rnode(9, mptrc) = hxposs(1)  ! dx
  rnode(10,mptrc) = hyposs(1)  ! dy
  rnode(11,mptrc) = possk(1)   ! dtbig
  rnode(4, mptrc) = FLOAT(node(6,mptrc)-1)*rnode(10,mptrc)
  rnode(5, mptrc) = FLOAT(node(5,mptrc)-1)*rnode(9,mptrc)
  rnode(6, mptrc) = rnode(4,mptrc)
  rnode(7, mptrc) = rnode(5,mptrc)
  rnode(8, mptrc) = rnode(2,mptrc)
  rnode(20,mptrc) = timestart
  rnode(29,mptrc) = 0.
  rnode(30,mptrc) = rnode(29,mptrc)                                     &
                  + FLOAT(node(14,mptrc)-1)*rnode(31,mptrc)

  CALL moment( rnode(1,mptrc),rnode(1,mptrc),4,usage,1.e+10,            &
               0.,rnode(9,mptrc),rnode(10,mptrc)            )

  RETURN
END SUBROUTINE domain
!

SUBROUTINE rdinit( a,b,m1,n1,k1,m,n,kk )
  DIMENSION a(m1,n1,k1),b(m,n,kk)
  DO ijk=1,m*n*kk
    b(ijk,1,1) = 0.
  END DO
  READ(12) a
  DO k=1,k1
    DO j=1,n1
      DO i=1,m1
        b(i,j,k) = a(i,j,k)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE rdinit
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

SUBROUTINE flcnst( mptr )
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
!
!  here we set the constants for the new grids
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)

!
! Store constant variables of the model into constant arrays
!
  CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)
!
  CALL resett
!
!  we're finished here
!
  RETURN
END SUBROUTINE flcnst
