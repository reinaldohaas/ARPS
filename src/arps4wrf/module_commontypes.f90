MODULE module_commontypes

  USE module_wrfgrid_constants

  TYPE TYPE_FIELD
!    integer :: ndims, istagger
!    integer, dimension(MAXDIMENSIONS) :: dom_start, mem_start, patch_start
!    integer, dimension(MAXDIMENSIONS) :: dom_end, mem_end, patch_end

    INTEGER  :: stagger
    INTEGER, DIMENSION(MAXDIMENSIONS) :: dom_start, mem_start, patch_start
    INTEGER, DIMENSION(MAXDIMENSIONS) :: dom_end,   mem_end,   patch_end

    REAL, POINTER, DIMENSION(:,:,:) :: rdata_arr

    CHARACTER (LEN=128), DIMENSION(MAXDIMENSIONS) :: dimnames
    CHARACTER (LEN=128) :: fieldname, mem_order, units, descr
  END TYPE TYPE_FIELD

END MODULE module_commontypes
