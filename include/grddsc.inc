!
!  the common block grddsc contains basic variable descriptions
!  the descriptions define:
!
!  stg*(*,*) -- grid staggering in the i'th coordinate direction
!            stg** -- ** denotes the 2-d or 3-d variable
!            (i,j) -- i denotes the i'th coordinate direction,
!                     and j denotes the variable number
!                     for the particular variable type
!                     xy, xz, yz, or xyz (see grid.h).
!                     staggering is in multiples of grid spacing.
!
!  idm*(*,*) -- storage similar to stg, the dimension of the
!            array for the variable will be n* - idm*(*,*).
!            we use arrays of all the same size, hence there is
!            some wasted space and this defines the true array
!            bounds
!
!  ipk*(*)   -- packing factor for storing grid data, storage
!            similar to stg*.
!
!  ini*(*,*) -- defines how to initialize arrays, (1,*) means
!            initialize, (2,*) is the other vector (see iup below)
!            (3,*) is the fields at dt earlier.
!
!  iup*(*,*) -- information on wether or not to update fields
!            in solution process.  if (1,*) is zero, no update.
!            if 1, update.  if (2,*) is nonzero, this is the
!            other component of this vector field and if it is
!            a negative number that means that this other
!            component is the first component (by righthand rule)
!
!  ibd*(*,*) -- information on wether or not to interpolate inner
!            grid boundary info to this field if necessary.
!            storage as in iup*(*,*)
!  ip*       -- contains pointers into the main storage array a
!            for the various grid fields and constants.
!            initialized in setstr.
!

  REAL :: stgxy(2,nxy2d  + 1),  stgxz(2,nxz2d  + 1),                    &
                    stgyz(2,nyz2d  + 1), stgxyz(3,nxyz3d + 1),          &
                    stgexbc(3,nexbc3d+ 1)

  INTEGER :: idmxy(2,nxy2d  + 1),  idmxz(2,nxz2d  + 1),                 &
                    idmyz(2,nyz2d  + 1), idmxyz(3,nxyz3d + 1),          &
                    idmexbc(3,nexbc3d+ 1),                              &
                    ipkxy(  nxy2d  + 1),  ipkxz(  nxz2d  + 1),          &
                    ipkyz(  nyz2d  + 1), ipkxyz(  nxyz3d + 1),          &
                    ipkexbc(nexbc3d+ 1),                                &
                    inixy(3,nxy2d  + 1),  inixz(3,nxz2d  + 1),          &
                    iniyz(3,nyz2d  + 1), inixyz(3,nxyz3d + 1),          &
                    iniexbc(3,nexbc3d+ 1),                              &
                    iupxy(2,nxy2d  + 1),  iupxz(2,nxz2d  + 1),          &
                    iupyz(2,nyz2d  + 1), iupxyz(2,nxyz3d + 1),          &
                    iupexbc(2,nexbc3d+ 1),                              &
                    ibdxy(3,nxy2d  + 1),  ibdxz(3,nxz2d  + 1),          &
                    ibdyz(3,nyz2d  + 1), ibdxyz(3,nxyz3d + 1),          &
                    ibdexbc(3,nexbc3d+ 1)

  COMMON / grddsc / stgxy,  stgxz, stgyz, stgxyz, stgexbc,              &
                    idmxy,  idmxz, idmyz, idmxyz, idmexbc,              &
                    ipkxy,  ipkxz, ipkyz, ipkxyz, ipkexbc,              &
                    inixy,  inixz, iniyz, inixyz, iniexbc,              &
                    iupxy,  iupxz, iupyz, iupxyz, iupexbc,              &
                    ibdxy,  ibdxz, ibdyz, ibdxyz, ibdexbc
!
!
!
!  this common area contains the pointers into storage for the
!  the grids
!
  INTEGER :: ipint(         30), ipreal(      30),  ips1d(         30), &
             ipx(           30), ipy(         30),  ipz(           30), &
             ipxy(nxy2d+1,  30), ipxz(nxz2d+1,30),  ipyz(nyz2d+1,  30), &
             ipxyz(nxyz3d+1,30), ipexbc(nexbc3d+1,30)

  COMMON / pntstr / ipint,ipreal,ips1d,                                 &
                    ipx,ipy,ipz,ipxy,ipxz,ipyz,ipxyz,ipexbc
!
!  next is a logical for a non-run storage check, (not used presently)
!
  LOGICAL :: checks
  COMMON /checkcm/ checks
