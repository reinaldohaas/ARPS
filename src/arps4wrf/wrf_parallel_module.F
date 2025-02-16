!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE PARALLEL_MODULE
!
! This module provides routines for parallelizing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module wrf_parallel_module

   integer, parameter :: IO_NODE = 0

   INTEGER :: HALO_WIDTH

   INTEGER :: ntprocs, my_proc_id, nxprocs, nyprocs

   integer :: my_x, my_y, &
              my_minx, my_miny, my_maxx, my_maxy

   integer, pointer, dimension(:,:) :: processors, &
                                       proc_minx, proc_maxx, &
                                       proc_miny, proc_maxy

   integer :: comm

#ifdef MPI_PARALLEL
   include 'mpif.h'
#endif

   contains


   SUBROUTINE wrf_parallel_init(myproc,nprocs_in,nprocx_in,nprocy_in)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: myproc
     INTEGER, INTENT(IN) :: nprocs_in, nprocx_in, nprocy_in

     INTEGER :: mpi_rank

     ntprocs = nprocs_in
     nxprocs = nprocx_in
     nyprocs = nprocy_in
     my_proc_id = myproc

#ifdef MPI_PARALLEL
     comm = MPI_COMM_WORLD
#else
     comm = 0
#endif

     mpi_rank = myproc

      ! Calculate which patch current processor will work on
     my_x = mod(mpi_rank,nxprocs)
     my_y = mpi_rank / nxprocs

     nullify(processors)
     nullify(proc_minx)
     nullify(proc_maxx)
     nullify(proc_miny)
     nullify(proc_maxy)

     RETURN
   end subroutine wrf_parallel_init

   SUBROUTINE wrf_parallel_set_halo_width(fake_width,istatus)

     IMPLICIT NONE
     INTEGER, INTENT(IN)  :: fake_width
     INTEGER, INTENT(OUT) :: istatus

     istatus = 0

     HALO_WIDTH = fake_width

     RETURN
   END SUBROUTINE wrf_parallel_set_halo_width

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_get_tile_dims
   !
   ! Purpose: To compute the starting and ending indices of the patch that the
   !   calling processor is to work on. When there are multiple processors,
   !   appropriate data structures describing the range of indices being
   !   worked on by other processors are also allocated and filled
   !   (processors, proc_minx, proc_maxx, proc_miny, proc_maxy).
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_get_tile_dims(xdim, ydim)

      implicit none

      ! Arguments
      integer, intent(in) :: xdim, ydim

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: i, j, ix, iy, px, py
      integer, dimension(2) :: buffer
      integer :: mpi_ierr
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      !
      ! Determine starting and ending grid points in x and y direction that we will work on
      !
      ! NOTE:
      ! For now, copy code from RSL_LITE's module_dm.F until build mechanism to link
      ! WRF and WPS code is worked out more.
      ! Eventually, it would probably be best to use module_dm code without copying
      !
      my_minx = -1
      j = 1
      do i = 1, xdim
         call task_for_point(i, j, 1, xdim, 1, ydim, nxprocs, nyprocs, px, py)
         if ( px == my_x ) then
            my_maxx = i
            if ( my_minx == -1 ) my_minx = i
         end if
      end do

      my_miny = -1
      i = 1
      do j = 1, ydim
         call task_for_point(i, j, 1, xdim, 1, ydim, nxprocs, nyprocs, px, py)
         if ( py == my_y ) then
            my_maxy = j
            if ( my_miny == -1 ) my_miny = j
         end if
      end do

      ! Create space to hold information about which other processors are
      !   working on which parts of the domain
      allocate(processors(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_minx(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_miny(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_maxx(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_maxy(0:nxprocs-1, 0:nyprocs-1))

      ! Exchange information with other processors
      if (my_proc_id == IO_NODE) then
         processors(my_x, my_y) = my_proc_id
         do i=1,ntprocs-1
            call MPI_Recv(buffer, 2, MPI_INTEGER, i, MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
            processors(buffer(1), buffer(2)) = mpi_stat(MPI_SOURCE)
         end do
      else
         buffer(1) = my_x
         buffer(2) = my_y
         call MPI_Send(buffer, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
      end if

      do ix=0,nxprocs-1
         do iy=0,nyprocs-1
            call parallel_bcast_int(processors(ix,iy), IO_NODE)
         end do
      end do

      proc_minx(my_x, my_y) = my_minx
      proc_maxx(my_x, my_y) = my_maxx
      proc_miny(my_x, my_y) = my_miny
      proc_maxy(my_x, my_y) = my_maxy

      do ix=0,nxprocs-1
         do iy=0,nyprocs-1
            call parallel_bcast_int(proc_minx(ix,iy), processors(ix,iy))
            call parallel_bcast_int(proc_maxx(ix,iy), processors(ix,iy))
            call parallel_bcast_int(proc_miny(ix,iy), processors(ix,iy))
            call parallel_bcast_int(proc_maxy(ix,iy), processors(ix,iy))
         end do
      end do

#else
      allocate(processors(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_minx(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_miny(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_maxx(0:nxprocs-1, 0:nyprocs-1))
      allocate(proc_maxy(0:nxprocs-1, 0:nyprocs-1))

      processors(0,0) = IO_NODE
      proc_minx(0,0) = 1
      proc_miny(0,0) = 1
      proc_maxx(0,0) = xdim
      proc_maxy(0,0) = ydim
      my_minx = 1
      my_maxx = xdim
      my_miny = 1
      my_maxy = ydim

#endif

   end subroutine parallel_get_tile_dims


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Copied from RSL_LITE's task_for_point.c until a good way can be found to
   !    get the build mechanism to use the original code in RSL_LITE.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine task_for_point(i_p, j_p, ids_p, ide_p, jds_p, jde_p, npx, npy, px, py)

      implicit none

      ! Arguments
      integer, intent(in) :: i_p, j_p, ids_p, ide_p, jds_p, jde_p, npx, npy
      integer, intent(out) :: px, py

      ! Local variables
      integer :: a, b, rem, idim, jdim, i, j, ids, jds, ide, jde

      i = i_p - 1
      j = j_p - 1
      ids = ids_p - 1
      jds = jds_p - 1
      ide = ide_p - 1
      jde = jde_p - 1

      idim = ide-ids+1
      jdim = jde-jds+1

      i = max(i,ids)
      i = min(i,ide)
      rem = mod(idim, npx)
      a = ( rem / 2 ) * ( (idim / npx) + 1 )
      b = a + ( npx - rem ) * ( idim / npx )
      if ( i-ids < a ) then
         px = (i-ids) / ( (idim / npx) + 1 )
      else if ( i-ids < b ) then
         px = ( a / ( (idim / npx) + 1 ) ) + (i-a-ids) / ( ( b - a ) / ( npx - rem ) )
      else
         px = ( a / ( (idim / npx) + 1 ) ) + (b-a-ids) / ( ( b - a ) / ( npx - rem ) ) + &
                                             (i-b-ids) / ( ( idim / npx ) + 1 )
      end if

      j = max(j,jds)
      j = min(j,jde)
      rem = mod(jdim, npy)
      a = ( rem / 2 ) * ( (jdim / npy) + 1 )
      b = a + ( npy - rem ) * ( jdim / npy )
      if ( j-jds < a ) then
         py = (j-jds) / ( (jdim / npy) + 1 )
      else if ( j-jds < b ) then
         py = ( a / ( (jdim / npy) + 1 ) ) + (j-a-jds) / ( ( b - a ) / ( npy - rem ) )
      else
         py = ( a / ( (jdim / npy) + 1 ) ) + (b-a-jds) / ( ( b - a ) / ( npy - rem ) ) + &
                                             (j-b-jds) / ( ( jdim / npy ) + 1 )
      end if

   end subroutine task_for_point


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: gather_whole_field_i
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gather_whole_field_i(patch_array, ms1, me1, ms2, me2, ms3, me3, &
                                 ps1, pe1, ps2, pe2, ps3, pe3,                &
                                 domain_array, ds1, de1, ds2, de2, ds3, de3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3, &
                             ds1, de1, ds2, de2, ds3, de3
      integer, dimension(ms1:me1,ms2:me2,ms3:me3), intent(in) :: patch_array
      integer, dimension(ds1:de1,ds2:de2,ds3:de3), intent(inout) :: domain_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: i, ii, j, jj, kk
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr, istatus
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
      INTEGER, ALLOCATABLE :: tem1(:,:,:)
      INTEGER :: nsize

      if (my_proc_id == IO_NODE) then

         do j=0,nyprocs-1
            do i=0,nxprocs-1
               if (processors(i,j) /= IO_NODE) then
                  call MPI_Recv(jdims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  call MPI_Recv(idims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  IF (ALLOCATED(tem1)) DEALLOCATE(tem1)
                  ALLOCATE(tem1( idims(1):idims(2),jdims(1):jdims(2), ds3:de3), STAT = istatus)
                  nsize = (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de3-ds3+1)
                  CALL MPI_Recv(tem1, nsize, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

                  do kk=ds3,de3
! BUG: Check on mpi_stat and mpi_ierr
!                     call MPI_Recv(domain_array(idims(1):idims(2),jdims(1):jdims(2),kk), &
!                                               (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1), &
!                                   MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                    DO jj = jdims(1), jdims(2)
                      DO ii = idims(1), idims(2)
                        domain_array(ii,jj,kk) = tem1(ii,jj,kk)
                      END DO
                    END DO
                  end do
               else
                  domain_array(ps1:pe1,ps2:pe2,ps3:pe3) = patch_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
         end do

      else

         jdims(1) = ps2
         jdims(2) = pe2
         call MPI_Send(jdims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call MPI_Send(idims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         ALLOCATE(tem1(ps1:pe1,ps2:pe2,ps3:pe3), STAT = istatus)
         nsize = (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1)
         do kk=ps3,pe3
!            call MPI_Send(patch_array(ps1:pe1,ps2:pe2,kk), (pe1-ps1+1)*(pe2-ps2+1), MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
! BUG: Check on mpi_ierr
          DO jj = ps2,pe2
            DO ii = ps1,pe1
              tem1(ii,jj,kk) = patch_array(ii,jj,kk)
            END DO
          END DO
         end do
         call MPI_Send(tem1, nsize, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
     end if

     DEALLOCATE(tem1)
#endif

   end subroutine gather_whole_field_i


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: gather_whole_field_r
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gather_whole_field_r(patch_array, ms1, me1, ms2, me2, ms3, me3, &
                                 ps1, pe1, ps2, pe2, ps3, pe3, &
                                 domain_array, ds1, de1, ds2, de2, ds3, de3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3, &
                             ds1, de1, ds2, de2, ds3, de3
      real, dimension(ms1:me1,ms2:me2,ms3:me3), intent(in) :: patch_array
      real, dimension(ds1:de1,ds2:de2,ds3:de3), intent(inout) :: domain_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: i, ii, j, jj, kk
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr, istatus
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      REAL, ALLOCATABLE :: tem1(:,:,:)
      INTEGER :: nsize

      if (my_proc_id == IO_NODE) then

         do j=0,nyprocs-1
            do i=0,nxprocs-1
               if (processors(i,j) /= IO_NODE) then
                  call MPI_Recv(jdims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  call MPI_Recv(idims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

                  IF (ALLOCATED(tem1)) DEALLOCATE(tem1)
                  ALLOCATE(tem1(idims(1):idims(2),jdims(1):jdims(2),ds3:de3), STAT = istatus)
                  nsize = (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de3-ds3+1)

                  CALL MPI_Recv(tem1,nsize, MPI_REAL, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

                  do kk=ds3,de3
! BUG: Check on mpi_stat and mpi_ierr
!                     call MPI_Recv(domain_array(idims(1):idims(2),jdims(1):jdims(2),kk), &
!                                   (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1), &
!                                   MPI_REAL, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                    DO jj = jdims(1),jdims(2)
                      DO ii = idims(1),idims(2)
                        domain_array(ii,jj,kk) = tem1(ii,jj,kk)
                      END DO
                    END DO
                  end do
               else
                  domain_array(ps1:pe1,ps2:pe2,ps3:pe3) = patch_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
         end do

      else

         jdims(1) = ps2
         jdims(2) = pe2
         call MPI_Send(jdims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call MPI_Send(idims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         ALLOCATE(tem1(idims(1):idims(2),jdims(1):jdims(2),ps3:pe3), STAT = istatus)
         nsize = (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(pe3-ps3+1)
         do kk=ps3,pe3
!            call MPI_Send(patch_array(ps1:pe1,ps2:pe2,kk), (pe1-ps1+1)*(pe2-ps2+1), MPI_REAL, 0, my_proc_id, comm, mpi_ierr)
! BUG: Check on mpi_ierr
           DO jj = ps2,pe2
            DO ii = ps1,pe1
              tem1(ii,jj,kk) = patch_array(ii,jj,kk)
            END DO
          END DO
         end do
         CALL MPI_SEND(tem1, nsize, MPI_REAL, 0, my_proc_id, comm, mpi_ierr)
      end if

      DEALLOCATE(tem1)
#endif

   end subroutine gather_whole_field_r


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: scatter_whole_field_i
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine scatter_whole_field_i(patch_array, ms1, me1, ms2, me2, ms3, me3, &
                                 ps1, pe1, ps2, pe2, ps3, pe3, &
                                 domain_array, ds1, de1, ds2, de2, ds3, de3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3, &
                             ds1, de1, ds2, de2, ds3, de3
      integer, dimension(ms1:me1,ms2:me2,ms3:me3), intent(inout) :: patch_array
      integer, dimension(ds1:de1,ds2:de2,ds3:de3), intent(in) :: domain_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: i, ii, j, jj, kk
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr, istatus
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      INTEGER, ALLOCATABLE :: tem1(:,:,:)
      INTEGER :: nsize

      if (my_proc_id == IO_NODE) then

         do j=0,nyprocs-1
            do i=0,nxprocs-1
               if (processors(i,j) /= IO_NODE) then
                  call MPI_Recv(jdims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  call MPI_Recv(idims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  IF (ALLOCATED(tem1)) DEALLOCATE(tem1)
                  ALLOCATE(tem1(idims(1):idims(2),jdims(1):jdims(2),ds3:de3), STAT = istatus)
                  nsize = (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de3-ds3+1)

                  do kk=ds3,de3
! BUG: Check on mpi_stat and mpi_ierr
!                     call MPI_Send(domain_array(idims(1):idims(2),jdims(1):jdims(2),kk), &
!                                               (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1), &
!                                   MPI_INTEGER, processors(i,j), my_proc_id, comm, mpi_ierr)
                    DO jj = jdims(1),jdims(2)
                      DO ii = idims(1),idims(2)
                        tem1(ii,jj,kk) = domain_array(ii,jj,kk)
                      END DO
                    END DO
                  end do

                  call MPI_Send(tem1,nsize, MPI_INTEGER, processors(i,j), my_proc_id, comm, mpi_ierr)

               else
                  patch_array(ps1:pe1,ps2:pe2,ps3:pe3) = domain_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
         end do

      else

         jdims(1) = ps2
         jdims(2) = pe2
         call MPI_Send(jdims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call MPI_Send(idims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)

         ALLOCATE(tem1(ps1:pe1,ps2:pe2,ps3:pe3), STAT = istatus)
         nsize = (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1)

         call MPI_Recv(tem1, nsize, MPI_INTEGER, 0, MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
!            call MPI_Recv(patch_array(ps1:pe1,ps2:pe2,kk), (pe1-ps1+1)*(pe2-ps2+1), &
!                          MPI_INTEGER, 0, MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
! BUG: Check on mpi_ierr
            DO jj = ps2,pe2
              DO ii = ps1,pe1
                patch_array(ii,jj,kk) = tem1(ii,jj,kk)
              END DO
            END DO

         end do
     end if

     DEALLOCATE(tem1)
#endif

   end subroutine scatter_whole_field_i


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: scatter_whole_field_r
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine scatter_whole_field_r(patch_array, ms1, me1, ms2, me2, ms3, me3, &
                                 ps1, pe1, ps2, pe2, ps3, pe3, &
                                 domain_array, ds1, de1, ds2, de2, ds3, de3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3, &
                             ds1, de1, ds2, de2, ds3, de3
      real, dimension(ms1:me1,ms2:me2,ms3:me3), intent(inout) :: patch_array
      real, dimension(ds1:de1,ds2:de2,ds3:de3), intent(in) :: domain_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: i, ii, j, jj, kk
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr, istatus
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      REAL, ALLOCATABLE :: tem1(:,:,:)
      INTEGER :: nsize

      if (my_proc_id == IO_NODE) then

         do j=0,nyprocs-1
            do i=0,nxprocs-1
               if (processors(i,j) /= IO_NODE) then
                  call MPI_Recv(jdims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
                  call MPI_Recv(idims, 2, MPI_INTEGER, processors(i,j), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

                  IF (ALLOCATED(tem1)) DEALLOCATE(tem1)
                  ALLOCATE(tem1(idims(1):idims(2),jdims(1):jdims(2),ds3:de3), STAT = istatus)
                  nsize = (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de3-ds3+1)

                  do kk=ds3,de3
! BUG: Check on mpi_stat and mpi_ierr
!                     call MPI_Send(domain_array(idims(1):idims(2),jdims(1):jdims(2),kk), &
!                                               (idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1), &
!                                   MPI_REAL, processors(i,j), my_proc_id, comm, mpi_ierr)
                    DO jj = jdims(1),jdims(2)
                      DO ii = idims(1),idims(2)
                        tem1(ii,jj,kk) = domain_array(ii,jj,kk)
                      END DO
                    END DO
                  end do

                  call MPI_Send(tem1,nsize, MPI_REAL, processors(i,j), my_proc_id, comm, mpi_ierr)
               else
                  patch_array(ps1:pe1,ps2:pe2,ps3:pe3) = domain_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
         end do

      else

         jdims(1) = ps2
         jdims(2) = pe2
         call MPI_Send(jdims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call MPI_Send(idims, 2, MPI_INTEGER, 0, my_proc_id, comm, mpi_ierr)

         ALLOCATE(tem1(ps1:pe1,ps2:pe2,ps3:pe3), STAT = istatus)
         nsize = (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1)

         call MPI_Recv(tem1, nsize, MPI_REAL, 0, MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
!            call MPI_Recv(patch_array(ps1:pe1,ps2:pe2,kk), (pe1-ps1+1)*(pe2-ps2+1), &
!                          MPI_REAL, 0, MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
! BUG: Check on mpi_ierr
            DO jj = ps2,pe2
              DO ii = ps1,pe1
                patch_array(ii,jj,kk) = tem1(ii,jj,kk)
              END DO
            END DO

         end do
     end if

     DEALLOCATE(tem1)
#endif

   end subroutine scatter_whole_field_r


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: exchange_halo_r
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine exchange_halo_r(patch_array, &
                              ms1, me1, ms2, me2, ms3, me3, &
                              ps1, pe1, ps2, pe2, ps3, pe3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3
      real, dimension(ms1:me1,ms2:me2,ms3:me3), intent(inout) :: patch_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: ii, jj, kk, i, j
      integer :: mpi_ierr
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      REAL, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:)
      INTEGER :: nsize
      INTEGER :: istatus

      IF (HALO_WIDTH <= 0) RETURN

      !
      ! Get left edge of halo
      !
      ALLOCATE(tem1(HALO_WIDTH,ms2:me2,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(HALO_WIDTH,ms2:me2,ps3:pe3), STAT = istatus)
      nsize = HALO_WIDTH*(me2-ms2+1)*(pe3-ps3+1)

      if (my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x+1,my_y), my_proc_id, comm, mpi_ierr)
              DO ii = pe1-HALO_WIDTH+1,pe1
                i = ii - pe1+HALO_WIDTH
                tem1(i,jj,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL,processors(my_x+1,my_y), my_proc_id, comm, mpi_ierr)
      end if

      if (my_x /= 0) then
         CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x-1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x-1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              DO ii = ms1,ms1+HALO_WIDTH-1
                i = ii - ms1 + 1
                patch_array(ii,jj,kk) = tem2(i,jj,kk)
              END DO
            end do
         end do
      end if

      !
      ! Get right edge of halo
      !
      if (my_x /= 0) then
         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x-1,my_y), my_proc_id, comm, mpi_ierr)
              DO ii = ps1,ps1+HALO_WIDTH-1
                i = ii - ps1 + 1
                tem1(i,jj,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL,processors(my_x-1,my_y), my_proc_id, comm, mpi_ierr)
      end if

      if (my_x /= (nxprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x+1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x+1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              DO ii = me1-HALO_WIDTH+1,me1
                i = ii - me1 + HALO_WIDTH
                patch_array(ii,jj,kk) = tem2(i,jj,kk)
              END DO
            end do
         end do
      end if

      DEALLOCATE(tem1, tem2)

      ALLOCATE(tem1(ms1:me1,HALO_WIDTH,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(ms1:me1,HALO_WIDTH,ps3:pe3), STAT = istatus)
      nsize = (me1-ms1+1)*HALO_WIDTH*(pe3-ps3+1)

      !
      ! Get bottom edge of halo
      !
      if (my_y /= (nyprocs - 1)) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_REAL, &
!                             processors(my_x,my_y+1), my_proc_id, comm, mpi_ierr)
              j = jj - pe2 + HALO_WIDTH
              DO ii = ms1,me1
                tem1(ii,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL, processors(my_x,my_y+1), my_proc_id, comm, mpi_ierr)
      end if

      if (my_y /= 0) then
         CALL MPI_RECV(tem2,nsize,MPI_REAL, processors(my_x,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!               call MPI_Recv(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_REAL, &
!                             processors(my_x,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              j = jj -ms2 + 1
              DO ii = ms1, me1
                patch_array(ii,jj,kk) = tem2(ii,j,kk)
              END DO
            end do
         end do
      end if

      !
      ! Get top edge of halo
      !
      if (my_y /= 0) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_REAL, &
!                             processors(my_x,my_y-1), my_proc_id, comm, mpi_ierr)
                j = jj - ps2 + 1
                DO ii = ms1,me1
                  tem1(ii,j,kk) = patch_array(ii,jj,kk)
                END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL,processors(my_x,my_y-1), my_proc_id, comm, mpi_ierr)
      end if

      if (my_y /= (nyprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!               call MPI_Recv(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_REAL, &
!                             processors(my_x,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = ms1, me1
                 patch_array(ii,jj,kk) = tem2(ii,j,kk)
               END DO
            end do
         end do
      end if
      DEALLOCATE(tem1, tem2)


      ALLOCATE(tem1(HALO_WIDTH,HALO_WIDTH,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(HALO_WIDTH,HALO_WIDTH,ps3:pe3), STAT = istatus)
      nsize = HALO_WIDTH*HALO_WIDTH*(pe3-ps3+1)

      !
      ! Get lower-right corner of halo
      !
      if (my_y /= (nyprocs - 1) .and. my_x /= 0) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x-1,my_y+1), my_proc_id, comm, mpi_ierr)
              j = jj - pe2 + HALO_WIDTH
              DO ii = ps1, ps1+HALO_WIDTH-1
                i = ii - ps1 + 1
                tem1(i,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL, processors(my_x-1,my_y+1), my_proc_id, comm, mpi_ierr)
      end if

      if (my_y /= 0 .and. my_x /= (nxprocs - 1)) then

        CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x+1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!               call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x+1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj -ms2 + 1
               DO ii = me1-HALO_WIDTH+1, me1
                 i = ii - me1 + HALO_WIDTH
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
            end do
         end do
      end if

      !
      ! Get upper-left corner of halo
      !
      if (my_y /= 0 .and. my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x+1,my_y-1), my_proc_id, comm, mpi_ierr)
              j = jj - ps2 + 1
              DO ii = pe1-HALO_WIDTH+1, pe1
                i = ii - pe1 + HALO_WIDTH
                tem1(i,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_REAL, processors(my_x+1,my_y-1), my_proc_id, comm, mpi_ierr)

      end if

      if (my_y /= (nyprocs - 1) .and. my_x /= 0) then

        CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x-1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!               call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x-1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = ms1, ms1+HALO_WIDTH-1
                 i = ii - ms1 + 1
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
              END DO
            end do
         end do
      end if

      !
      ! Get upper-right corner of halo
      !
      if (my_y /= 0 .and. my_x /= 0) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x-1,my_y-1), my_proc_id, comm, mpi_ierr)
               j = jj -ps2+1
               DO ii = ps1, ps1+HALO_WIDTH-1
                 i = ii - ps1 +1
                 tem1(i,j,kk) = patch_array(ii,jj,kk)
               END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_REAL,processors(my_x-1,my_y-1), my_proc_id, comm, mpi_ierr)
      end if

      if (my_y /= (nyprocs - 1) .and. my_x /= (nxprocs - 1)) then

        CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x+1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!              call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                            processors(my_x+1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = me1-HALO_WIDTH+1, me1
                 i = ii - me1 + HALO_WIDTH
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
            end do
         end do
      end if

      !
      ! Get lower-left corner of halo
      !
      if (my_y /= (nyprocs - 1) .and. my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                             processors(my_x+1,my_y+1), my_proc_id, comm, mpi_ierr)
                j = jj - pe2+ HALO_WIDTH
                DO ii = pe1-HALO_WIDTH+1,pe1
                  i = ii - pe1 + HALO_WIDTH
                  tem1(i,j,kk) = patch_array(ii,jj,kk)
                END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_REAL, processors(my_x+1,my_y+1), my_proc_id, comm, mpi_ierr)
      end if

      if (my_y /= 0 .and. my_x /= 0) then
        CALL MPI_RECV(tem2,nsize,MPI_REAL,processors(my_x-1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!              call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_REAL, &
!                            processors(my_x-1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - ms2 + 1
               DO ii = ms1,ms1+HALO_WIDTH-1
                 i = ii - ms1 + 1
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
            end do
         end do
      end if

      DEALLOCATE(tem1, tem2)
#endif

   end subroutine exchange_halo_r


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: exchange_halo_i
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine exchange_halo_i(patch_array, &
                              ms1, me1, ms2, me2, ms3, me3, &
                              ps1, pe1, ps2, pe2, ps3, pe3)

      implicit none

      ! Arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ms1, me1, ms2, me2, ms3, me3
      integer, dimension(ms1:me1,ms2:me2,ms3:me3), intent(inout) :: patch_array

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: ii, jj, kk, i, j
      integer :: mpi_ierr, istatus
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
      INTEGER, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:)
      INTEGER :: nsize

      IF (HALO_WIDTH <= 0) RETURN

      ALLOCATE(tem1(HALO_WIDTH,ms2:me2,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(HALO_WIDTH,ms2:me2,ps3:pe3), STAT = istatus)
      nsize = HALO_WIDTH*(me2-ms2+1)*(pe3-ps3+1)

      !
      ! Get left edge of halo
      !
      if (my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y), my_proc_id, comm, mpi_ierr)
              DO ii = pe1-HALO_WIDTH+1,pe1
                i = ii - pe1+HALO_WIDTH
                tem1(i,jj,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do

         CALL MPI_SEND(tem1,nsize,MPI_INTEGER,processors(my_x+1,my_y), my_proc_id, comm, mpi_ierr)
      end if

      if (my_x /= 0) then

         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x-1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              DO ii = ms1,ms1+HALO_WIDTH-1
                i = ii - ms1 + 1
                patch_array(ii,jj,kk) = tem2(i,jj,kk)
              END DO

            end do
         end do
      end if

      !
      ! Get right edge of halo
      !
      if (my_x /= 0) then
         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y), my_proc_id, comm, mpi_ierr)
              DO ii = ps1,ps1+HALO_WIDTH-1
                i = ii - ps1 + 1
                tem1(i,jj,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER,processors(my_x-1,my_y), my_proc_id, comm, mpi_ierr)

      end if
      if (my_x /= (nxprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x+1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)

         do kk=ps3,pe3
            do jj=ms2,me2
!               call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              DO ii = me1-HALO_WIDTH+1,me1
                i = ii - me1 + HALO_WIDTH
                patch_array(ii,jj,kk) = tem2(i,jj,kk)
              END DO
            end do
         end do
      end if
      DEALLOCATE(tem1, tem2)

      ALLOCATE(tem1(ms1:me1,HALO_WIDTH,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(ms1:me1,HALO_WIDTH,ps3:pe3), STAT = istatus)
      nsize = (me1-ms1+1)*HALO_WIDTH*(pe3-ps3+1)

      !
      ! Get bottom edge of halo
      !
      if (my_y /= (nyprocs - 1)) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_INTEGER, &
!                             processors(my_x,my_y+1), my_proc_id, comm, mpi_ierr)
              j = jj - pe2 + HALO_WIDTH
              DO ii = ms1,me1
                tem1(ii,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER, processors(my_x,my_y+1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= 0) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER, processors(my_x,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!               call MPI_Recv(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_INTEGER, &
!                             processors(my_x,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
              j = jj -ms2 + 1
              DO ii = ms1, me1
                patch_array(ii,jj,kk) = tem2(ii,j,kk)
              END DO
            end do
         end do
      end if

      !
      ! Get top edge of halo
      !
      if (my_y /= 0) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_INTEGER, &
!                             processors(my_x,my_y-1), my_proc_id, comm, mpi_ierr)
                j = jj - ps2 + 1
                DO ii = ms1,me1
                  tem1(ii,j,kk) = patch_array(ii,jj,kk)
                END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER,processors(my_x,my_y-1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= (nyprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!               call MPI_Recv(patch_array(ms1:me1,jj,kk), (me1-ms1+1), MPI_INTEGER, &
!                             processors(my_x,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = ms1, me1
                 patch_array(ii,jj,kk) = tem2(ii,j,kk)
               END DO
            end do
         end do
      end if
      DEALLOCATE(tem1, tem2)


      ALLOCATE(tem1(HALO_WIDTH,HALO_WIDTH,ps3:pe3), STAT = istatus)
      ALLOCATE(tem2(HALO_WIDTH,HALO_WIDTH,ps3:pe3), STAT = istatus)
      nsize = HALO_WIDTH*HALO_WIDTH*(pe3-ps3+1)
      !
      ! Get lower-right corner of halo
      !
      if (my_y /= (nyprocs - 1) .and. my_x /= 0) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y+1), my_proc_id, comm, mpi_ierr)
              j = jj - pe2 + HALO_WIDTH
              DO ii = ps1, ps1+HALO_WIDTH-1
                i = ii - ps1 + 1
                tem1(i,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER, processors(my_x-1,my_y+1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= 0 .and. my_x /= (nxprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x+1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!               call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj -ms2 + 1
               DO ii = me1-HALO_WIDTH+1, me1
                 i = ii - me1 + HALO_WIDTH
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
            end do
         end do
      end if

      !
      ! Get upper-left corner of halo
      !
      if (my_y /= 0 .and. my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y-1), my_proc_id, comm, mpi_ierr)
              j = jj - ps2 + 1
              DO ii = pe1-HALO_WIDTH+1, pe1
                i = ii - pe1 + HALO_WIDTH
                tem1(i,j,kk) = patch_array(ii,jj,kk)
              END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER, processors(my_x+1,my_y-1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= (nyprocs - 1) .and. my_x /= 0) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x-1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!               call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = ms1, ms1+HALO_WIDTH-1
                 i = ii - ms1 + 1
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
              END DO
            end do
         end do
      end if

      !
      ! Get upper-right corner of halo
      !
      if (my_y /= 0 .and. my_x /= 0) then
         do kk=ps3,pe3
            do jj=ps2,ps2+HALO_WIDTH-1
!               call MPI_Send(patch_array(ps1:ps1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y-1), my_proc_id, comm, mpi_ierr)
               j = jj -ps2+1
               DO ii = ps1, ps1+HALO_WIDTH-1
                 i = ii - ps1 +1
                 tem1(i,j,kk) = patch_array(ii,jj,kk)
               END DO
            end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER,processors(my_x-1,my_y-1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= (nyprocs - 1) .and. my_x /= (nxprocs - 1)) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x+1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=me2-HALO_WIDTH+1,me2
!               call MPI_Recv(patch_array(me1-HALO_WIDTH+1:me1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y+1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - me2 + HALO_WIDTH
               DO ii = me1-HALO_WIDTH+1, me1
                 i = ii - me1 + HALO_WIDTH
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
             end do
         end do
      end if

      !
      ! Get lower-left corner of halo
      !
      if (my_y /= (nyprocs - 1) .and. my_x /= (nxprocs - 1)) then
         do kk=ps3,pe3
            do jj=pe2-HALO_WIDTH+1,pe2
!               call MPI_Send(patch_array(pe1-HALO_WIDTH+1:pe1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x+1,my_y+1), my_proc_id, comm, mpi_ierr)
                j = jj - pe2+ HALO_WIDTH
                DO ii = pe1-HALO_WIDTH+1,pe1
                  i = ii - pe1 + HALO_WIDTH
                  tem1(i,j,kk) = patch_array(ii,jj,kk)
                END DO
              end do
         end do
         CALL MPI_SEND(tem1,nsize,MPI_INTEGER, processors(my_x+1,my_y+1), my_proc_id, comm, mpi_ierr)
      end if
      if (my_y /= 0 .and. my_x /= 0) then
         CALL MPI_RECV(tem2,nsize,MPI_INTEGER,processors(my_x-1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
         do kk=ps3,pe3
            do jj=ms2,ms2+HALO_WIDTH-1
!               call MPI_Recv(patch_array(ms1:ms1+HALO_WIDTH-1,jj,kk), HALO_WIDTH, MPI_INTEGER, &
!                             processors(my_x-1,my_y-1), MPI_ANY_TAG, comm, mpi_stat, mpi_ierr)
               j = jj - ms2 + 1
               DO ii = ms1,ms1+HALO_WIDTH-1
                 i = ii - ms1 + 1
                 patch_array(ii,jj,kk) = tem2(i,j,kk)
               END DO
            end do
         end do
      end if

      DEALLOCATE(tem1, tem2)
#endif

   end subroutine exchange_halo_i


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_bcast_logical
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_bcast_logical(lval)

      implicit none

      ! Argument
      logical, intent(inout) :: lval

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      call MPI_Bcast(lval, 1, MPI_LOGICAL, IO_NODE, comm, mpi_ierr)
#endif

   end subroutine parallel_bcast_logical


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_bcast_int
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_bcast_int(ival, from_whom)

      implicit none

      ! Argument
      integer, intent(inout) :: ival
      integer, intent(in), optional :: from_whom

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      if (present(from_whom)) then
         call MPI_Bcast(ival, 1, MPI_INTEGER, from_whom, comm, mpi_ierr)
      else
         call MPI_Bcast(ival, 1, MPI_INTEGER, IO_NODE, comm, mpi_ierr)
      end if
#endif

   end subroutine parallel_bcast_int


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_bcast_real
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_bcast_real(rval, from_whom)

      implicit none

      ! Argument
      real, intent(inout) :: rval
      integer, intent(in), optional :: from_whom

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      if (present(from_whom)) then
         call MPI_Bcast(rval, 1, MPI_REAL, from_whom, comm, mpi_ierr)
      else
         call MPI_Bcast(rval, 1, MPI_REAL, IO_NODE, comm, mpi_ierr)
      end if
#endif

   end subroutine parallel_bcast_real


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_bcast_char
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_bcast_char(cval, n, from_whom)

      implicit none

      ! Argument
      integer, intent(in) :: n
      character (len=n), intent(inout) :: cval
      integer, intent(in), optional :: from_whom

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      if (present(from_whom)) then
         call MPI_Bcast(cval, n, MPI_CHARACTER, from_whom, comm, mpi_ierr)
      else
         call MPI_Bcast(cval, n, MPI_CHARACTER, IO_NODE, comm, mpi_ierr)
      end if
#endif

   end subroutine parallel_bcast_char


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_finish
   !
   ! Purpose: Free up, deallocate, and for MPI, finalize.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_finish()

      implicit none

      ! Arguments

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      call MPI_Finalize(mpi_ierr)
#endif

      if (associated(processors)) deallocate(processors)
      if (associated(proc_minx)) deallocate(proc_minx)
      if (associated(proc_maxx)) deallocate(proc_maxx)
      if (associated(proc_miny)) deallocate(proc_miny)
      if (associated(proc_maxy)) deallocate(proc_maxy)

   end subroutine parallel_finish


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_abort
   !
   ! Purpose: Terminate everything
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_abort()

      implicit none

      ! Arguments

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr, mpi_errcode

      call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#endif

      stop

   end subroutine parallel_abort

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: parallel_bcast_real_array
   !
   ! Purpose: Broadcast a real array (new subroutine from CAPS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parallel_bcast_real_array(rval,nsize, from_whom)

      implicit none

      ! Argument
      integer, intent(in) :: nsize
      real, intent(inout) :: rval(nsize)
      integer, intent(in), optional :: from_whom

      ! Local variables
#ifdef MPI_PARALLEL
      integer :: mpi_ierr

      if (present(from_whom)) then
         call MPI_Bcast(rval, nsize, MPI_REAL, from_whom, comm, mpi_ierr)
      else
         call MPI_Bcast(rval, nsize, MPI_REAL, IO_NODE, comm, mpi_ierr)
      end if
#endif

   end subroutine parallel_bcast_real_array
end module wrf_parallel_module
