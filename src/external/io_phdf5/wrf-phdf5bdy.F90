!* The HDF5 WRF IO module was written by the the HDF Group at NCSA, the     *
!* National Center for Supercomputing Applications.                         *
!*     HDF Group                                                            *
!*     National Center for Supercomputing Applications                      *
!*     University of Illinois at Urbana-Champaign                           *
!*     605 E. Springfield, Champaign IL 61820                               *
!*     http://hdf.ncsa.uiuc.edu/                                            *
!*                                                                          *
!* Copyright 2004 by the Board of Trustees, University of Illinois,         *
!*                                                                          *
!* Redistribution or use of this IO module, with or without modification,   *
!* is permitted for any purpose, including commercial  purposes.            *
!*                                                                          *
!* This software is an unsupported prototype.  Use at your own risk.        *
!*     http://hdf.ncsa.uiuc.edu/apps/WRF-ROMS                               *
!*                                                                          *
!* This work was funded by the MEAD expedition at the National Center       *
!* for Supercomputing Applications, NCSA.  For more information see:        *
!*     http://www.ncsa.uiuc.edu/expeditions/MEAD                            *
!*                                                                          *
!*                                                                          *
!****************************************************************************/
!
!##################################################################
!######                                                      ######
!######     This file was modified specifically for          ######
!######                 WRF boundary file                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
!  Changes:
!
!     1. HDF5IOWRITE_bdy write dataset independently based on a flag
!        - IOFlag, which is passed in from ext_phdf5_write_bdy.
!     2. ext_phdf5_write_bdy is the same as ext_phdf5_write_field except
!        that it passes an extra LOGICAL dummy parameter, IOFlag.
!     3. Set no_par = .FALSE. i.e. does not support compact data storage
!       
!-----------------------------------------------------------------------

subroutine HDF5IOWRITE_bdy(DataHandle,Comm,DateStr,Length,DomainStart,DomainEnd &
     ,PatchStart,PatchEnd,MemoryOrder &
     ,WrfDType,FieldType,groupID,TimeIndex &
     ,DimRank ,DatasetName,XField,IOFlag,Status)

  use wrf_phdf5_data
  use ext_phdf5_support_routines
  use HDF5
  implicit none
  include 'mpif.h'
  include 'wrf_status_codes.h'

  integer                     ,intent(in)     :: DataHandle
  integer                     ,intent(inout)  :: Comm
  character*(*)               ,intent(in)     :: DateStr
  integer,dimension(NVarDims) ,intent(in)     :: Length

  integer,dimension(NVarDims) ,intent(in)     :: DomainStart
  integer,dimension(NVarDims) ,intent(in)     :: DomainEnd
  integer,dimension(NVarDims) ,intent(in)     :: PatchStart
  integer,dimension(NVarDims) ,intent(in)     :: PatchEnd

  character*(*)               ,intent(in)     :: MemoryOrder

  integer                     ,intent(in)     :: WrfDType
  integer(hid_t)              ,intent(in)     :: FieldType
  integer(hid_t)              ,intent(in)     :: groupID
  integer                     ,intent(in)     :: TimeIndex

  integer,dimension(*)        ,intent(in)     :: DimRank
  character (*)               ,intent(in)     :: DatasetName
  integer,dimension(*)        ,intent(inout)  :: XField
  LOGICAL,                     INTENT(IN)     :: IOFlag
  integer                     ,intent(out)    :: Status

  integer(hid_t)                              :: dset_id
  integer                                     :: NDim
  integer,dimension(NVarDims)                 :: VStart
  integer,dimension(NVarDims)                 :: VCount
  character (3)                               :: Mem0
  character (3)                               :: UCMem0
  type(wrf_phdf5_data_handle) ,pointer         :: DH

  ! attribute defination
  integer(hid_t)                              :: dimaspace_id  ! DimRank dataspace id
  integer(hid_t)                              :: dimattr_id    ! DimRank attribute id
  integer(hsize_t) ,dimension(1)              :: dim_space
  INTEGER(HID_T)                              :: dspace_id     ! Raw Data memory Dataspace id
  INTEGER(HID_T)                              :: fspace_id     ! Raw Data file Dataspace id
  INTEGER(HID_T)                              :: crp_list      ! chunk  identifier
  integer(hid_t)                              :: h5_atypeid    ! for fieldtype,memorder attribute
  integer(hid_t)                              :: h5_aspaceid   ! for fieldtype,memorder  
  integer(hid_t)                              :: h5_attrid     ! for fieldtype,memorder
  integer(hsize_t), dimension(7)              :: adata_dims
  integer                                     :: routine_atype


  integer,          dimension(:),allocatable  :: dimrank_data

  INTEGER(HSIZE_T), dimension(:),allocatable  :: dims  ! Dataset dimensions
  INTEGER(HSIZE_T), dimension(:),allocatable  :: sizes ! Dataset dimensions
  INTEGER(HSIZE_T), dimension(:),allocatable  :: chunk_dims 
  INTEGER(HSIZE_T), dimension(:),allocatable  :: hdf5_maxdims
  INTEGER(HSIZE_T), dimension(:),allocatable  :: offset 
  INTEGER(HSIZE_T), dimension(:),allocatable  :: count  
  INTEGER(HSIZE_T), DIMENSION(7)              :: dimsfi
  integer                                     :: hdf5err
  integer                                     :: i,j
  integer(size_t)                             :: dsetsize

  ! FOR PARALLEL IO
  integer(hid_t)                              :: xfer_list
  logical                                     :: no_par

!-------- uncomment for debug ----------------
!
!  INTEGER :: myproc
!
!  CALL mpi_comm_rank( Comm, myproc, Status )
!
!-------- uncomment for debug ----------------

  ! get the handle 
  call GetDH(DataHandle,DH,Status)
  if(Status /= WRF_NO_ERR) then
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     return
  endif

  ! get the rank of the dimension
  call GetDim(MemoryOrder,NDim,Status)
  if(Status /= WRF_NO_ERR) then
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     return
  endif

  ! If patch is equal to domain, the parallel is not necessary, sequential is used.
  ! In this version, we haven't implemented this yet.
  ! We use no_par to check whether we can use compact data storage.
!
! CAPS begin:
!
!  no_par = .TRUE.
!  do i = 1,NDim
!     if((PatchStart(i)/=DomainStart(i)).or.(PatchEnd(i)/=DomainEnd(i))) then
!        no_par = .FALSE.
!        exit
!     endif
!  enddo
!
!-----------------------------------------------------------------------
!
!  This is not a correct criteria. anyway, we do not use H5D_COMPACT 
!  data storage for bounary data.
!
!-----------------------------------------------------------------------

  no_par = .FALSE.
!
! CAPS end
!

  ! change the different Memory Order to XYZ for patch and domain
  if(MemoryOrder.NE.'0') then
     call ExtOrder(MemoryOrder,PatchStart,Status)
     call ExtOrder(MemoryOrder,PatchEnd,Status)
     call ExtOrder(MemoryOrder,DomainStart,Status)
     call ExtOrder(MemoryOrder,DomainEnd,Status)
  endif

  ! allocating memory for dynamic arrays; 
  ! since the time step is always 1, we may ignore the fourth
  ! dimension time; now keep it to make it consistent with sequential version
  allocate(dims(NDim+1))
  allocate(count(NDim+1))
  allocate(offset(NDim+1))
  allocate(sizes(NDim+1))


  ! arrange offset, count for each hyperslab
  dims(1:NDim)   = DomainEnd(1:NDim) - DomainStart(1:NDim) + 1
  dims(NDim+1)   = 1

  count(NDim+1)  = 1
  count(1:NDim)  = Length(1:NDim)

  offset(NDim+1) = 0
  offset(1:NDim) = PatchStart(1:NDim) - 1

  ! allocate the dataspace to write hyperslab data

  dimsfi = 0
  do i = 1, NDim + 1
     dimsfi(i) = count(i)
  enddo

  ! create the memory space id
  call h5screate_simple_f(NDim+1,count,dspace_id,hdf5err,count)
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_DATASPACE
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif


  ! create file space
  call h5screate_simple_f(NDim+1,dims,fspace_id,hdf5err,dims)
  if(hdf5err.lt.0) then        
     Status =  WRF_HDF5_ERR_DATASPACE
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  ! compact storage when the patch is equal to the whole domain
  ! calculate the non-decomposed dataset size

  call h5tget_size_f(FieldType,dsetsize,hdf5err)
  if(hdf5err.lt.0) then
     Status = WRF_HDF5_ERR_DATATYPE
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  do i =1,NDim
     dsetsize = dsetsize*dims(i)
  enddo
  if(no_par.and.(dsetsize.le.CompDsetSize)) then
     call h5pcreate_f(H5P_DATASET_CREATE_F,crp_list,hdf5err)
     if(hdf5err.lt.0) then
        Status =  WRF_HDF5_ERR_PROPERTY_LIST
        write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
        call wrf_debug ( WARN , msg) 
        deallocate(dims)
        deallocate(count)
        deallocate(offset)
        deallocate(sizes)
        return
     endif
     call h5pset_layout_f(crp_list,0,hdf5err)
     if(hdf5err.lt.0) then
        Status =  WRF_HDF5_ERR_PROPERTY_LIST
        write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
        call wrf_debug ( WARN , msg) 
        deallocate(dims)
        deallocate(count)
        deallocate(offset)
        deallocate(sizes)
        return
     endif
     call h5dcreate_f(DH%TgroupIDs(TimeIndex),DatasetName,FieldType,fspace_id,dset_id,&
          hdf5err,crp_list)
     call h5pclose_f(crp_list,hdf5err)
  else
     call h5dcreate_f(DH%TgroupIDs(TimeIndex),DatasetName,FieldType,fspace_id,dset_id,hdf5err)
  endif

  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_DATASET_CREATE 
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  ! select the correct hyperslab for file space id
  CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, count &
       ,hdf5err) 
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_DATASET_GENERAL
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  ! Create property list for collective dataset write
  CALL h5pcreate_f(H5P_DATASET_XFER_F, xfer_list, hdf5err)
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_PROPERTY_LIST
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  CALL h5pset_dxpl_mpio_f(xfer_list, H5FD_MPIO_COLLECTIVE_F&
       ,hdf5err)
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_PROPERTY_LIST
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  ! write the data in memory space to file space
!-----------------------------------------------------------------------
!
! CAPS Begin: 
!   Write dataset independently instead of colletively
!
!-----------------------------------------------------------------------

!  CALL h5dwrite_f(dset_id,FieldType,XField,dimsfi,hdf5err,&
!       mem_space_id =dspace_id,file_space_id =fspace_id, &
!       xfer_prp = xfer_list)

  IF (IOFlag) THEN

!write(0,*) ':: ',myproc,'  is writing ',datasetname

    CALL h5dwrite_f(dset_id,FieldType,XField,dimsfi,hdf5err,&
       mem_space_id =dspace_id,file_space_id =fspace_id)
    if(hdf5err.lt.0) then
       Status =  WRF_HDF5_ERR_DATASET_WRITE
       write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
       call wrf_debug ( WARN , msg) 
       deallocate(dims)
       deallocate(count)
       deallocate(offset)
       deallocate(sizes)
       return
    endif

  END IF
!-----------------------------------------------------------------------
!
! CAPS End: 
!
!-----------------------------------------------------------------------

  CALL h5pclose_f(xfer_list,hdf5err)
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_PROPERTY_LIST
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif

  if(TimeIndex == 1) then
     do i =1, MaxVars
        if(DH%dsetids(i) == -1) then
           DH%dsetids(i) = dset_id
           DH%VarNames(i) = DataSetName
           exit
        endif
     enddo
     ! Only writing attributes when TimeIndex ==1
     call write_hdf5_attributes(DataHandle,MemoryOrder,WrfDType,DimRank,&
          NDim,dset_id,Status)
  endif

  call h5sclose_f(fspace_id,hdf5err)
  call h5sclose_f(dspace_id,hdf5err)
  if(TimeIndex /= 1) then
     call h5dclose_f(dset_id,hdf5err)  
  endif
  if(hdf5err.lt.0) then
     Status =  WRF_HDF5_ERR_DATASPACE  
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg) 
     deallocate(dims)
     deallocate(count)
     deallocate(offset)
     deallocate(sizes)
     return
  endif
  Status = WRF_NO_ERR
  return
end subroutine  HDF5IOWRITE_bdy

! The real routine to write HDF5 boundary file
subroutine ext_phdf5_write_bdy(DataHandle,DateStr,Var,Field,FieldType,&
     Comm,IOComm,DomainDesc,MemoryOrder,  &
     Stagger,DimNames,DomainStart,DomainEnd,&
     MemoryStart,MemoryEnd,PatchStart,PatchEnd,&
     IOFlag,Status)

  use wrf_phdf5_data
  use ext_phdf5_support_routines
  USE HDF5 ! This module contains all necessary modules 
  implicit none
  include 'wrf_status_codes.h'

  integer                       ,intent(in)      :: DataHandle
  character*(*)                 ,intent(in)      :: DateStr
  character*(*)                 ,intent(in)      :: Var
  integer                       ,intent(inout)   :: Field(*)
  integer                       ,intent(in)      :: FieldType
  integer                       ,intent(inout)   :: Comm
  integer                       ,intent(inout)   :: IOComm
  integer                       ,intent(in)      :: DomainDesc
  character*(*)                 ,intent(in)      :: MemoryOrder
  character*(*)                 ,intent(in)      :: Stagger   ! Dummy for now
  character*(*) , dimension (*) ,intent(in)      :: DimNames
  integer ,dimension(*)         ,intent(in)      :: DomainStart, DomainEnd
  integer ,dimension(*)         ,intent(in)      :: MemoryStart, MemoryEnd
  integer ,dimension(*)         ,intent(in)      :: PatchStart,  PatchEnd
  LOGICAL,                       INTENT(IN)      :: IOFlag
  integer                       ,intent(out)     :: Status

  type(wrf_phdf5_data_handle)    ,pointer        :: DH
  integer(hid_t)                                 :: GroupID
  integer                                        :: NDim
  character (VarNameLen)                         :: VarName
  character (3)                                  :: MemO
  character (3)                                  :: UCMemO
  integer(hid_t)                                 :: DsetID
  integer      ,dimension(NVarDims)              :: Length
  integer      ,dimension(NVarDims)              :: DomLength
  integer      ,dimension(NVarDims+1)            :: DimRank
  character(256),dimension(NVarDims)              :: RODimNames
  integer      ,dimension(NVarDims)              :: StoredStart
  integer      ,dimension(:,:,:,:),allocatable   :: XField
  integer      ,dimension(:,:,:,:),allocatable   :: BUFFER! for logical field
  integer                                        :: stat
  integer                                        :: NVar
  integer                                        :: i,j,k,m,dim_flag
  integer                                        :: i1,i2,j1,j2,k1,k2
  integer                                        :: x1,x2,y1,y2,z1,z2
  integer                                        :: l1,l2,m1,m2,n1,n2
  integer(hid_t)                                 :: XType
  integer                                        :: di
  character (256)                                 :: NullName
  integer                                        :: TimeIndex
  integer ,dimension(NVarDims+1)                 :: temprank
  logical                                        :: NotFound


  NullName = char(0)
  dim_flag = 0

  call GetDH(DataHandle,DH,Status)
  if(Status /= WRF_NO_ERR) then
     write(msg,*) 'Warning Status = ',Status,' in ',__FILE__,', line', __LINE__
     call wrf_debug ( WARN , msg)
     return
  endif

  ! Examine here, Nov. 7th, 2003
  if(DH%FileStatus == WRF_FILE_OPENED_AND_COMMITTED) then 

     ! obtain group id and initialize the rank of dimensional attributes
     GroupID = DH%GroupID
     DimRank = -1

     ! get the rank of the dimension based on MemoryOrder string(cleaver from NetCDF)
     call GetDim(MemoryOrder,NDim,Status)
     if(Status /= WRF_NO_ERR) then
        write(msg,*) 'Warning BAD MEMORY ORDER in ',__FILE__,', line', __LINE__
        call wrf_debug ( WARN , msg)
        return
     endif

     ! check whether the DateStr is the correct length
     call DateCheck(DateStr,Status)
     if(Status /= WRF_NO_ERR) then
        write(msg,*) 'Warning DATE STRING ERROR in ',__FILE__,', line', __LINE__ 
        call wrf_debug ( WARN , msg)
        return
     endif

     ! get the dataset name and dimensional information of the data
     VarName           = Var
     Length(1:NDim)    = PatchEnd(1:NDim) - PatchStart(1:NDim) + 1
     DomLength(1:NDim) = DomainEnd(1:NDim) - DomainStart(1:NDim) + 1

     ! Transposing the data order and dim. string order, store to RODimNames
     call ExtOrder(MemoryOrder,Length,Status)
     call ExtOrder(MemoryOrder,DomLength,Status)
     if(Status /= WRF_NO_ERR) then
        write(msg,*) 'Warning BAD MEMORY ORDER in ',__FILE__,', line', __LINE__ 
        call wrf_debug ( WARN , msg)
        return
     endif

     ! Map datatype from WRF to HDF5
     select case (FieldType)
     case (WRF_REAL)
        XType = H5T_NATIVE_REAL
     case (WRF_REAL8)
        Xtype = H5T_NATIVE_DOUBLE
     case (WRF_INTEGER)
        XType = H5T_NATIVE_INTEGER
     case (WRF_LOGICAL)
        XType = DH%EnumID
     case default
        Status = WRF_HDF5_ERR_DATA_TYPE_NOTFOUND
        return
     end select

     ! HANDLE  with dim. scale 
     ! handle dimensional scale data; search and store them in a table.
     ! The table is one dimensional array of compound data type. One member of
     ! the type is HDF5 string, representing the name of the dim(west_east_stag eg.)
     ! Another number is the length of the dimension(west_east_stag = 31)
     ! In this part, we will not store TIME but leave it at the end since the time
     ! index won't be known until the end of the run; since all fields(HDF5 datasets)
     ! have the same timestamp, writing it once should be fine.

     ! 1) create a loop for dimensions
     call GetDataTimeIndex('write',DataHandle,DateStr,TimeIndex,Status)
     if(Status /= WRF_NO_ERR) then
        return
     endif

     if(TimeIndex == 1) then

        ! 2) get the dim. name, the first dim. is reserved for time, 
        call ExtOrderStr(MemoryOrder,DimNames,RODimNames,Status)
        if(Status /= WRF_NO_ERR) then
           write(msg,*) 'Warning BAD MEMORY ORDER in ',__FILE__,', line', __LINE__ 
           call wrf_debug ( WARN , msg)
           return
        endif
        ! 3) get the dim. length
        ! 4) inside the loop, search the table for dimensional name( table module)
        !    IF FOUND, go to the next dimension, return the table dimensional rank
        !    (For example, find west_east_stag in the table, the rank of "west_east_stag"
        !     is 3; so return 3 for the array dimrank.)
        !    in the table; so through the table, we can find the information
        !    such as names, length of this dimension
        ! 4.1) save the rank into an array for attribute
        !      if not found,  go to 5)
        ! 4)' the first dimension is reserved for time, so table starts from j = 2 
        !
        ! 5) NOT FOUND, inside the loop add the new dimensional information to the 
        ! table(table module)

        ! The first dimension of the field is always "time" and "time"
        ! is also the first dimension of the "table".
        k = 2
        DimRank(1) = 1

        do i = 1,NDim
           do j = 2,MaxTabDims

              ! Search for the table and see if we are at the end of the table
              if (DH%DIMTABLE(j)%dim_name == NO_NAME) then

                 ! Sometimes the RODimNames is NULLName or ''. If that happens,
                 ! we will search the table from the beginning and see 
                 ! whether the name is FAKEDIM(the default name) and  the 
                 ! current length of the dim. is the same as that of FAKEDIM; 
                 ! if yes, use this FAKEDIM for the current field dim. 

                 if(RODimNames(i) ==''.or. RODimNames(i)==NullName) then
                    do m = 2,j
                       if(DomLength(i)==DH%DIMTABLE(m)%Length.and. &
                            DH%DIMTABLE(m)%dim_name(1:7)=='FAKEDIM')then
                          DimRank(k) = m
                          k = k + 1
                          dim_flag = 1
                          exit
                       endif
                    enddo
                    ! No FAKEDIM and the same length dim. is found,
                    ! Add another dimension "FAKEDIM + j", with the length
                    ! as DomLength(i)
                    if (dim_flag == 1) then 
                       dim_flag = 0
                    else   
                       RODimNames(i) = 'FAKEDIM'//achar(j+iachar('0'))
                       DH%DIMTABLE(j)%dim_name  = RODimNames(i)
                       DH%DIMTABLE(j)%length    = DomLength(i)
                       DimRank(k) = j
                       k          = k + 1
                    endif
                    ! no '' or NULLName is found, then assign this RODimNames
                    ! to the dim. table.
                 else
                    DH%DIMTABLE(j)%dim_name  = RODimNames(i)
                    DH%DIMTABLE(j)%length    = DomLength(i)
                    DimRank(k)               = j
                    k = k + 1
                 endif
                 exit
                 ! If we found the current dim. in the table already,save the rank
              else if(DH%DIMTABLE(j)%dim_name == RODimNames(i)) then
                 ! remember the rank of dimensional scale
                 DimRank(k) = j
                 k = k + 1
                 exit
              else
                 continue
              endif
           enddo
        enddo
     endif ! end of timeindex of 1

     ! 6) create an attribute array called DimRank to store the rank of the attribute.
     !    This will be done in the HDF5IOWRITE routine        

     ! 7) before the end of the run, 1) update time, 2) write the table to HDF5.

     ! get the index of l1,.......for writing HDF5 file.
     StoredStart = 1
     call GetIndices(NDim,MemoryStart,MemoryEnd,l1,l2,m1,m2,n1,n2)
     call GetIndices(NDim,StoredStart,Length   ,x1,x2,y1,y2,z1,z2)
     call GetIndices(NDim,PatchStart, PatchEnd ,i1,i2,j1,j2,k1,k2)
     di=1
     if(FieldType == WRF_REAL8) di = 2
     allocate(XField(di,x1:x2,y1:y2,z1:z2), STAT=stat)
     if(stat/= 0) then
        Status = WRF_ERR_FATAL_ALLOCATION_ERROR
        write(msg,*) 'Fatal ALLOCATION ERROR in ',__FILE__,', line',__LINE__
        call wrf_debug ( FATAL , msg)
        return
     endif

     ! Transpose the real data for tools people
     call Transpose_hdf5('write',MemoryOrder,di, Field,l1,l2,m1,m2,n1,n2 &
          ,XField,x1,x2,y1,y2,z1,z2 &
          ,i1,i2,j1,j2,k1,k2 )

     ! handle with logical data separately,because of not able to 
     ! map Fortran Logical type to C type
     if(FieldType .eq. WRF_LOGICAL) then
        allocate(BUFFER(di,x1:x2,y1:y2,z1:z2), STAT=stat)
        do k =z1,z2
           do j = y1,y2
              do i = x1,x2
                 do m = 1,di
                    if(XField(m,i,j,k)/= 0) then
                       BUFFER(m,i,j,k) = 1
                    else
                       BUFFER(m,i,j,k) = 0
                    endif
                 enddo
              enddo
           enddo
        enddo
        call HDF5IOWRITE_bdy(DataHandle,Comm,DateStr,Length,DomainStart, DomainEnd &
             ,PatchStart,PatchEnd, MemoryOrder &
             ,FieldType,XType,groupID,TimeIndex,DimRank &
             ,Var,BUFFER,IOFlag,Status)
        deallocate(BUFFER,STAT=stat)
        if(stat/=0) then
           Status = WRF_ERR_FATAL_ALLOCATION_ERROR
           write(msg,*) 'Fatal ALLOCATION ERROR in ',__FILE__,', line',__LINE__
           call wrf_debug ( FATAL , msg)
           return
        endif
     else 
        call HDF5IOWRITE_bdy(DataHandle,Comm,DateStr,Length, DomainStart, DomainEnd &
             ,PatchStart, PatchEnd, MemoryOrder &
             ,FieldType,XType,groupID,TimeIndex,DimRank &
             ,Var,XField,IOFlag,Status)
     endif

     if (Status /= WRF_NO_ERR) then 
        return
     endif

     deallocate(XField,STAT=stat)
     if(stat/=0) then
        Status = WRF_ERR_FATAL_ALLOCATION_ERROR
        write(msg,*) 'Fatal ALLOCATION ERROR in ',__FILE__,', line',__LINE__
        call wrf_debug ( FATAL , msg)
        return
     endif
  endif

  return 

end subroutine ext_phdf5_write_bdy
