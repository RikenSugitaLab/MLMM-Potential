!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_pairlist_str_mod
!> @brief   structure of pairlist information
!! @authors Yuji Sugita (YS), Takashi Imai (TI), Jaewoon Jung (JJ), 
!!          Takaharu Mori (TM), Motoshi Kamiya (MK), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_pairlist_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_table_pair
    integer                       :: num_nb15_max
    integer                       :: water_nb15_max
    integer,          allocatable :: num_nb15_calc(:)
    integer,          allocatable :: num_nb15_calcw(:)
    integer,          allocatable :: nb15_calc_list(:,:)
    integer,          allocatable :: nb15_calc_listw(:,:)
    integer,          allocatable :: cell_linked_list(:)
    integer,          allocatable :: cell_linked_listw(:)
    integer,          allocatable :: atom_cell_index(:)
    integer,          allocatable :: num_nb15_calc_water(:)
    integer,          allocatable :: nb15_calc_list_water(:,:)
    integer,          allocatable :: cell_linked_list_water(:)
    integer,          allocatable :: atom_cell_index_water(:)
    integer,          allocatable :: num_list_water(:)
  end type s_table_pair

  type, public :: s_pairlist
    type(s_table_pair)            :: table
    logical                       :: allocation
    logical                       :: allocate_nobc
    logical                       :: allocate_pbc
    logical                       :: allocate_solsol
    logical                       :: allocate_solwat
    logical                       :: allocate_watwat
    integer                       :: num_nb15_max
    integer,          allocatable :: num_nb15_pre(:)
    integer,          allocatable :: num_nb15(:)
    integer,          allocatable :: num_nb15_calc(:,:)
    integer,          allocatable :: nb15_calc_list(:,:)
    real(wp)                      :: pairlistdist
    integer                       :: ecqm_num_nb15
    integer, allocatable          :: ecqm_nb15_list(:,:)
    ! for GBSA
    integer                       :: num_all_max
    integer,          allocatable :: num_all_pre(:)
    integer,          allocatable :: num_all(:)
    integer,          allocatable :: num_all_calc(:,:)
    integer,          allocatable :: all_calc_list(:,:)
  end type s_pairlist

  ! parameters for allocatable variables
  integer,      public, parameter :: PairListAtomNobc        = 1
  integer,      public, parameter :: PairListIntNobc         = 2
  integer,      public, parameter :: PairListPbcSolute       = 3
  integer,      public, parameter :: PairListPbcWater        = 4
  integer,      public, parameter :: PairListPbcSoluteSolute = 5
  integer,      public, parameter :: PairListPbcSoluteWater  = 6
  integer,      public, parameter :: PairListPbcWaterWater   = 7
  integer,      public, parameter :: PairListNthreads        = 8
  integer,      public, parameter :: PairListEcqm            = 9
  integer,      public, parameter :: PairListAtomNobcGbsa    = 10
  integer,      public, parameter :: PairListIntNobcGbsa     = 11
  integer,      public, parameter :: PairListNthreadsGbsa    = 12

  ! subroutines
  public :: init_pairlist
  public :: alloc_pairlist
  public :: alloc_pairlist2
  public :: dealloc_pairlist
  public :: dealloc_pairlist_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_pairlist
  !> @brief        initialize pairlist information
  !! @authors      YS
  !! @param[out]   pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_pairlist(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    pairlist%table%num_nb15_max   = 0
    pairlist%table%water_nb15_max = 0

    pairlist%allocation           = .true.
    pairlist%allocate_nobc        = .true.
    pairlist%allocate_pbc         = .true.
    pairlist%allocate_solsol      = .true.
    pairlist%allocate_solwat      = .true.
    pairlist%allocate_watwat      = .true.
    pairlist%num_nb15_max         = 0
    pairlist%pairlistdist         = 0.0_wp

    return

  end subroutine init_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist
  !> @brief        allocate pairlist information
  !! @authors      YS, TI, JJ, TM
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist(pairlist, variable, var_size)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: nthread, omp_get_max_threads


    alloc_stat   = 0
    dealloc_stat = 0

    ! get number of threads
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! allocate selected variables
    !
    select case (variable)

    case (PairListAtomNobc)

      if (allocated(pairlist%num_nb15_calc)) then
        if (size(pairlist%num_nb15_calc) == var_size*nthread) return
        deallocate(pairlist%num_nb15_calc, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%num_nb15_calc(var_size,nthread), &
               stat = alloc_stat)

    case (PairListIntNobc)

      if (allocated(pairlist%nb15_calc_list)) then
        if (size(pairlist%nb15_calc_list) == var_size*nthread) return
        deallocate(pairlist%nb15_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%nb15_calc_list(var_size,nthread), &
               stat = alloc_stat)

    case (PairListPbcSolute)

      if (allocated(pairlist%table%num_nb15_calc)) then
        if (size(pairlist%table%num_nb15_calc) == var_size) return
        deallocate(pairlist%table%num_nb15_calc,    &
                   pairlist%table%num_nb15_calcw,   &
                   pairlist%table%atom_cell_index,  &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%num_nb15_calc(var_size),    &
               pairlist%table%num_nb15_calcw(var_size),   &
               pairlist%table%cell_linked_list(var_size), &
               pairlist%table%atom_cell_index(var_size),  &
               stat = alloc_stat)

    case (PairListPbcWater)

      if (allocated(pairlist%table%num_nb15_calc_water)) then
        if (size(pairlist%table%num_nb15_calc_water) == var_size) return
        deallocate(pairlist%table%num_nb15_calc_water,    &
                   pairlist%table%cell_linked_list_water, &
                   pairlist%table%atom_cell_index_water,  &
                   pairlist%table%cell_linked_listw,      &
                   pairlist%table%num_list_water,         &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%num_nb15_calc_water(var_size),    &
               pairlist%table%cell_linked_list_water(var_size), &
               pairlist%table%atom_cell_index_water(var_size),  &
               pairlist%table%cell_linked_listw(3*var_size),    &
               pairlist%table%num_list_water(var_size),         &
               stat = alloc_stat)

    case (PairListNthreads)

      if (allocated(pairlist%num_nb15_pre)) then
        if (size(pairlist%num_nb15_pre) == nthread) return
        deallocate(pairlist%num_nb15_pre, &
                   pairlist%num_nb15,     &
                   stat = dealloc_stat)
      end if
      allocate(pairlist%num_nb15_pre(nthread), &
               pairlist%num_nb15    (nthread), &
               stat = alloc_stat)

    case (PairListEcqm)

      if (allocated(pairlist%ecqm_nb15_list)) then
        if (size(pairlist%ecqm_nb15_list) == var_size*2) return
        deallocate(pairlist%ecqm_nb15_list, stat = dealloc_stat)
      end if
      allocate(pairlist%ecqm_nb15_list(2, var_size), stat = alloc_stat)

    case (PairListAtomNobcGbsa)

      if (allocated(pairlist%num_all_calc)) then
        if (size(pairlist%num_all_calc) == var_size*nthread) return
        deallocate(pairlist%num_all_calc, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%num_all_calc(var_size,nthread), &
               stat = alloc_stat)

    case (PairListIntNobcGbsa)

      if (allocated(pairlist%all_calc_list)) then
        if (size(pairlist%all_calc_list) == var_size*nthread) return
        deallocate(pairlist%all_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%all_calc_list(var_size,nthread), &
               stat = alloc_stat)

    case (PairListNthreadsGbsa)

      if (allocated(pairlist%num_all_pre)) then
        if (size(pairlist%num_all_pre) == nthread) return
        deallocate(pairlist%num_all_pre, &
                   pairlist%num_all,     &
                   stat = dealloc_stat)
      end if
      allocate(pairlist%num_all_pre(nthread), &
               pairlist%num_all    (nthread), &
               stat = alloc_stat)

    case  default

      call error_msg('Alloc_PairList> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist2
  !> @brief        allocate pairlist information
  !! @authors      MK
  !! @param[inout] pairlist  : pairlist information
  !! @param[in]    variable  : allocatable variables
  !! @param[in]    var_size1 : size of variables (1st dimension)
  !! @param[in]    var_size2 : size of variables (2nd dimension)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist2(pairlist, variable, var_size1, var_size2)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: nthread, omp_get_max_threads


    alloc_stat   = 0
    dealloc_stat = 0

    ! get number of threads
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! allocate selected variables
    !
    select case (variable)

    case (PairListPbcSoluteSolute)

      if (allocated(pairlist%table%nb15_calc_list)) then
        if (size(pairlist%table%nb15_calc_list) == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_list(var_size1,var_size2), &
               stat = alloc_stat)

    case (PairListPbcSoluteWater)

      if (allocated(pairlist%table%nb15_calc_listw)) then
        if (size(pairlist%table%nb15_calc_listw) == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_listw, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_listw(var_size1,var_size2), &
               stat = alloc_stat)

    case (PairListPbcWaterWater)

      if (allocated(pairlist%table%nb15_calc_list_water)) then
        if (size(pairlist%table%nb15_calc_list_water)               &
            == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_list_water, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_list_water(var_size1,var_size2), &
               stat = alloc_stat)

    case  default

      call error_msg('Alloc_PairList2> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_pairlist2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist
  !> @brief        deallocate pairlist information
  !! @authors      YS, TI, JJ
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist(pairlist, variable)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat

  
    dealloc_stat = 0

    select case (variable)

    case (PairListAtomNobc)

      if (allocated(pairlist%num_nb15_calc)) then
        deallocate (pairlist%num_nb15_calc,   &
                    stat = dealloc_stat)
      end if

    case (PairListIntNobc)

      if (allocated(pairlist%nb15_calc_list)) then
        deallocate (pairlist%nb15_calc_list, &
                    stat = dealloc_stat)
      end if

    case (PairListPbcSolute)

      if (allocated(pairlist%table%num_nb15_calc)) then
        deallocate (pairlist%table%num_nb15_calc,    &
                    pairlist%table%num_nb15_calcw,   &
                    pairlist%table%cell_linked_list, &
                    pairlist%table%atom_cell_index,  &
                    stat = dealloc_stat)
      end if

    case (PairListPbcWater)

      if (allocated(pairlist%table%num_nb15_calc_water)) then
        deallocate (pairlist%table%num_nb15_calc_water,    &
                    pairlist%table%cell_linked_list_water, &
                    pairlist%table%atom_cell_index_water,  &
                    stat = dealloc_stat)
      end if

    case (PairListPbcSoluteSolute)

      if (allocated(pairlist%table%nb15_calc_list)) then
        deallocate( pairlist%table%nb15_calc_list, &
                    stat = dealloc_stat)
      end if

    case (PairListPbcSoluteWater)

      if (allocated(pairlist%table%nb15_calc_listw)) then
        deallocate (pairlist%table%nb15_calc_listw, &
                    stat = dealloc_stat)
      end if

    case (PairListPbcWaterWater)

      if (allocated(pairlist%table%nb15_calc_list_water)) then
        deallocate (pairlist%table%nb15_calc_list_water, &
                    stat = dealloc_stat)
      end if

    case (PairListNthreads)

      if (allocated(pairlist%num_nb15_pre)) then
        deallocate (pairlist%num_nb15_pre,     &
                    pairlist%num_nb15,         &
                    stat = dealloc_stat)
      end if

    case (PairListEcqm)

      if (allocated(pairlist%ecqm_nb15_list)) then
        deallocate (pairlist%ecqm_nb15_list, stat = dealloc_stat)
      end if

    case (PairListAtomNobcGbsa)

      if (allocated(pairlist%num_all_calc)) then
        deallocate (pairlist%num_all_calc,   &
                    stat = dealloc_stat)
      end if

    case (PairListIntNobcGbsa)

      if (allocated(pairlist%all_calc_list)) then
        deallocate (pairlist%all_calc_list, &
                    stat = dealloc_stat)
      end if

    case (PairListNthreadsGbsa)

      if (allocated(pairlist%num_all_pre)) then
        deallocate (pairlist%num_all_pre,     &
                    pairlist%num_all,         &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_PairList> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist_all
  !> @brief        deallocate all pairlist information
  !! @authors      YS, TI
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist_all(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    call dealloc_pairlist(pairlist, PairListAtomNobc)
    call dealloc_pairlist(pairlist, PairListIntNobc)
    call dealloc_pairlist(pairlist, PairListPbcSolute)
    call dealloc_pairlist(pairlist, PairListPbcWater)
    call dealloc_pairlist(pairlist, PairListPbcSoluteSolute)
    call dealloc_pairlist(pairlist, PairListPbcSoluteWater)
    call dealloc_pairlist(pairlist, PairListPbcWaterWater)
    call dealloc_pairlist(pairlist, PairListNthreads)
    call dealloc_pairlist(pairlist, PairListEcqm)
    call dealloc_pairlist(pairlist, PairListAtomNobcGbsa)
    call dealloc_pairlist(pairlist, PairListIntNobcGbsa)
    call dealloc_pairlist(pairlist, PairListNthreadsGbsa)

    return

  end subroutine dealloc_pairlist_all

end module at_pairlist_str_mod
