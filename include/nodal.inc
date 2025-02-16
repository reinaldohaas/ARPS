!
  REAL :: rnode,tol,bzone,hxposs,hyposs,possk
  INTEGER :: node,lstart,newstl,llist,lback,mstart,ndfree,              &
          intrat,intratt,lfine,kcheck,mxnest,ncheck,levmlt

  COMMON / nodal / rnode(35,30),node(25,30),lstart(10),newstl(10),      &
                   llist(2,200),lback(10),tol,bzone,mstart,ndfree,      &
                   intrat,intratt,lfine,kcheck,mxnest,                  &
                   hxposs(10),hyposs(10),possk(10),                     &
                   ncheck(10),levmlt(10)
!
!   data definitions for node and rnode in nodal.h
!
!   grid information for grid mptr is stored in array rnode(*,mptr)
!   and node(*,mptr).  the information stored is
!
!   rnode(i,mptr), i =
!
!  1) x - corner 1
!  2) y - corner 1
!  3) x - corner 2
!  4) y - corner 2
!  5) x - corner 3
!  6) y - corner 3
!  7) x - corner 4
!  8) y - corner 4
!  9) dx
!    10) dy
!    11) dt
!    12) e vec1 x
!    13) e vec1 y
!    14) e vec2 x
!    15) e vec2 y
!    16) e max1
!    17) e min1
!    18) e max2
!    19) e min2
!    20) time to which grid has been integrated
!    21) cos of orientation angle
!      using x axis.  x axis corner 1 to corner 4.  coarse
!      grid (mptr=1) has cos=1,sin=0.
!    22) sin of orientation angle
!    23) slope corners 1 to 2
!    24) slope corners 1 to 4
!    25) x intercept corners 1 to 2
!    26) x intercept corners 1 to 4
!    27) x intercept corners 3 to 4
!    28) x intercept corners 2 to 3
!    29) grid bottom height (computational space - excluding terrain)
!    30) grid top height
!    31) dz
!
!   node(i,mptr), i =
!
!  1) parent grid
!  2) kid grid
!  3) sibling
!  4) nestlevel
!  5) number of rows (nx)
!  6) number of columns (ny)
!  7)
!  8)
!  9) pointer into intersect list (llist)
!    10) next grid on level pointer
!    11)
!    12) previous grid on level pointer
!    13) number of timesteps taken on grid
!    14) number of vertical levels (nz)
!    15)
