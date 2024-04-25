!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_str_mod
!> @brief   structure of pairlist information
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_str_mod

  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pairlist
    real(wp)                        :: pairlistdist
    integer,            allocatable :: num_nb15_calc(:,:)
    integer,            allocatable :: num_nb15_calcs(:,:)
    integer,            allocatable :: num_nb15_calcw(:,:)
    integer,            allocatable :: num_nb15_calcww(:,:)
    integer,            allocatable :: num_nb15_calc1(:,:)
    integer,            allocatable :: num_nb15_water(:,:)
    integer,            allocatable :: num_nb15_nobc(:,:)
    integer,            allocatable :: nb15_calc_list(:,:)
    integer,            allocatable :: nb15_calc_list1(:,:)
    integer,            allocatable :: nb15_calc_water(:,:)
    integer,            allocatable :: nb15_calc_list_nobc(:,:,:,:)
    integer,            allocatable :: nb15_cell(:)
    integer,            allocatable :: nb15_cells(:)
    integer,            allocatable :: nb15_cellw(:)
    integer,            allocatable :: nb15_cellww(:)
    integer,            allocatable :: nb15_list(:,:)
    integer,            allocatable :: nb15_lists(:,:)
    integer,            allocatable :: nb15_listw(:,:)
    integer,            allocatable :: nb15_listww(:,:)
    ! for GPU
    integer,            allocatable :: univ_ij_load(:)         !  (ij)
    ! for FEP
    integer,            allocatable :: num_nb15_calc_fep(:,:)
    integer,            allocatable :: num_nb15_calc1_fep(:,:)
    integer,            allocatable :: nb15_calc_list_fep(:,:)
    integer,            allocatable :: nb15_calc_list1_fep(:,:)
    integer,            allocatable :: nb15_cell_fep(:)
    integer,            allocatable :: nb15_list_fep(:,:)

#ifndef PGICUDA
    integer(1),         allocatable :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer(1),         allocatable :: pack_univ_mask2(:)      !  (ix*iy*ij/8)
    integer,            allocatable :: univ_ij_sort_list(:)    !  (ij)
    integer,            allocatable :: univ_ix_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_ix_list(:,:)       !  (iix, ij)
    integer,            allocatable :: univ_iy_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_iy_list(:,:)       !  (iiy, ij)
#else
    integer(1), allocatable, pinned :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer,    allocatable, pinned :: univ_ij_sort_list(:)    !  (ij)
    integer,    allocatable, pinned :: univ_ix_natom(:)        !  (ij)
    integer(1), allocatable, pinned :: univ_ix_list(:,:)       !  (iix, ij)
    integer,    allocatable, pinned :: univ_iy_natom(:)        !  (ij)
    integer(1), allocatable, pinned :: univ_iy_list(:,:)       !  (iiy, ij)
#endif
    integer                 :: univ_ncell_nonzero
    integer                 :: univ_update = 0
    integer                 :: univ_mask2_size = 0
    integer                 :: pack_univ_mask2_size = 0
  end type s_pairlist

  ! parameters for allocatable variables
  integer,        public, parameter :: PairListNoTable     = 1
  integer,        public, parameter :: PairListTable       = 2
  integer,        public, parameter :: PairListNoTable_CG  = 3
  integer,        public, parameter :: PairListNOBC        = 4
  integer,        public, parameter :: PairListNoTable_FEP = 5

  ! variables for maximum numbers in one cell
  integer,        public            :: MaxNb15      = 15000
  integer,        public            :: MaxNb15_NOBC = 3000
  integer,        public            :: MaxNb15_CG   = 3000
  integer,        public            :: MaxNb15_chk  = 15000
  integer,        public            :: MaxNb15Water = 1500
  ! FEP
  ! MaxNb15_fep will be changed depending on the numbe of 
  ! perturbed atoms.
  integer,        public            :: MaxNb15_fep  = 15000

  ! subroutines
  public  :: init_pairlist
  public  :: alloc_pairlist
  public  :: dealloc_pairlist
  public  :: dealloc_pairlist_all

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


    pairlist%pairlistdist = 0.0_wp

    return

  end subroutine init_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist
  !> @brief        allocate pairlist information
  !! @authors      JJ
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
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    MaxNb15_chk = MaxNB15

    ! allocate selected variables
    !
    select case (variable)

    case (PairListNOBC)

      if (allocated(pairlist%num_nb15_nobc)) then
        if (size(pairlist%num_nb15_nobc) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_nobc,       &
                     pairlist%nb15_calc_list_nobc, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_nobc))                              &
        allocate(pairlist%num_nb15_nobc(MaxAtom, var_size),                     &
                 pairlist%nb15_calc_list_nobc(2,MaxNb15_NOBC,MaxAtom,var_size), &
                 stat = alloc_stat)

      pairlist%num_nb15_nobc(1:MaxAtom, 1:var_size)   = 0

    case (PairListNoTable)

#ifndef USE_GPU
      if (allocated(pairlist%num_nb15_calc1)) then
        if (size(pairlist%num_nb15_calc1) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1,  &
                     pairlist%num_nb15_calc,   &
                     pairlist%nb15_calc_list1, &
                     pairlist%nb15_calc_list,  &
                     pairlist%nb15_cell,       &
                     pairlist%nb15_list,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1)) &
        allocate(pairlist%num_nb15_calc1 (MaxAtom, var_size), &
                 pairlist%num_nb15_calc  (MaxAtom, maxcell),  &
                 pairlist%nb15_calc_list1(MaxNb15, var_size), &
                 pairlist%nb15_calc_list (MaxNb15, maxcell),  &
                 pairlist%nb15_cell      (maxcell),           &
                 pairlist%nb15_list      (MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%num_nb15_calc1 (1:MaxAtom, 1:var_size)  = 0
      pairlist%num_nb15_calc  (1:MaxAtom, 1:maxcell)   = 0
#else
      if (.not. allocated(pairlist%num_nb15_calc1)) then
         ! debug
        allocate(pairlist%univ_ij_load(univ_maxcell1),            &
                 pairlist%univ_ij_sort_list(univ_maxcell1),       &
                 pairlist%univ_ix_natom(univ_maxcell1),           &
                 pairlist%univ_ix_list(MaxAtom, univ_maxcell1),   &
                 pairlist%univ_iy_natom(univ_maxcell1),           &
                 pairlist%univ_iy_list(MaxAtom, univ_maxcell1),   &
                 stat = alloc_stat)
      endif
      call set_pinned_memory(pairlist%univ_ij_sort_list, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_ix_natom, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_iy_natom, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_ix_list, MaxAtom*univ_maxcell1)
      call set_pinned_memory(pairlist%univ_iy_list, MaxAtom*univ_maxcell1)

      ! initialization
      pairlist%univ_ij_load(:)      = 0
      pairlist%univ_ij_sort_list(:) = 0
      pairlist%univ_ix_natom(:)     = 0
      pairlist%univ_ix_list(:,:)    = 0
      pairlist%univ_iy_natom(:)     = 0
      pairlist%univ_iy_list(:,:)    = 0
#endif


    case (PairListNoTable_FEP)

#ifndef USE_GPU
      if (allocated(pairlist%num_nb15_calc1)) then
        if (size(pairlist%num_nb15_calc1) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1,        &
                     pairlist%num_nb15_calc,         &
                     pairlist%nb15_calc_list1,       &
                     pairlist%nb15_calc_list,        &
                     pairlist%nb15_cell,             &
                     pairlist%nb15_list,             &
                     pairlist%num_nb15_calc1_fep,  &
                     pairlist%num_nb15_calc_fep,   &
                     pairlist%nb15_calc_list1_fep, &
                     pairlist%nb15_calc_list_fep,  &
                     pairlist%nb15_cell_fep,       &
                     pairlist%nb15_list_fep,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1)) &
        allocate(pairlist%num_nb15_calc1   (MaxAtom, var_size),     &
                 pairlist%num_nb15_calc    (MaxAtom, maxcell),      &
                 pairlist%nb15_calc_list1  (MaxNb15, var_size),     &
                 pairlist%nb15_calc_list   (MaxNb15, maxcell),      &
                 pairlist%nb15_cell        (maxcell),               &
                 pairlist%nb15_list        (MaxAtom, maxcell),      &
                 pairlist%num_nb15_calc1_fep (MaxAtom, var_size), &
                 pairlist%num_nb15_calc_fep  (MaxAtom, maxcell),  &
                 pairlist%nb15_calc_list1_fep(MaxNb15_fep, var_size), &
                 pairlist%nb15_calc_list_fep (MaxNb15_fep, maxcell),  &
                 pairlist%nb15_cell_fep      (maxcell),           &
                 pairlist%nb15_list_fep      (MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%num_nb15_calc1   (1:MaxAtom, 1:var_size)      = 0
      pairlist%num_nb15_calc    (1:MaxAtom, 1:maxcell)       = 0
      pairlist%num_nb15_calc1_fep (1:MaxAtom, 1:var_size)  = 0
      pairlist%num_nb15_calc_fep  (1:MaxAtom, 1:maxcell)   = 0
#else
      if (.not. allocated(pairlist%num_nb15_calc1)) then
         ! debug
        allocate(pairlist%univ_ij_load(univ_maxcell1),            &
                 pairlist%univ_ij_sort_list(univ_maxcell1),       &
                 pairlist%univ_ix_natom(univ_maxcell1),           &
                 pairlist%univ_ix_list(MaxAtom, univ_maxcell1),   &
                 pairlist%univ_iy_natom(univ_maxcell1),           &
                 pairlist%univ_iy_list(MaxAtom, univ_maxcell1),   &
                 stat = alloc_stat)
      endif
      call set_pinned_memory(pairlist%univ_ij_sort_list, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_ix_natom, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_iy_natom, univ_maxcell1*4)
      call set_pinned_memory(pairlist%univ_ix_list, MaxAtom*univ_maxcell1)
      call set_pinned_memory(pairlist%univ_iy_list, MaxAtom*univ_maxcell1)

      ! initialization
      pairlist%univ_ij_load(:)      = 0
      pairlist%univ_ij_sort_list(:) = 0
      pairlist%univ_ix_natom(:)     = 0
      pairlist%univ_ix_list(:,:)    = 0
      pairlist%univ_iy_natom(:)     = 0
      pairlist%univ_iy_list(:,:)    = 0

#endif


    case (PairListTable)

      if (allocated(pairlist%num_nb15_calc1)) then
        if (size(pairlist%num_nb15_calc1) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1,  &
                     pairlist%num_nb15_calc,   &
                     pairlist%num_nb15_calcs,  &
                     pairlist%num_nb15_calcw,  &
                     pairlist%num_nb15_calcww, &
                     pairlist%num_nb15_water,  &
                     pairlist%nb15_cell,       &
                     pairlist%nb15_cells,      &
                     pairlist%nb15_cellw,      &
                     pairlist%nb15_cellww,     &
                     pairlist%nb15_list,       &
                     pairlist%nb15_lists,      &
                     pairlist%nb15_listw,      &
                     pairlist%nb15_listww,     &
                     pairlist%nb15_calc_list1, &
                     pairlist%nb15_calc_list,  &
                     pairlist%nb15_calc_water, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1)) &
        allocate(pairlist%num_nb15_calc1 (MaxAtom, var_size),     &
                 pairlist%num_nb15_calc  (MaxAtom, maxcell),      &
                 pairlist%num_nb15_calcs (MaxAtom, maxcell),      &
                 pairlist%num_nb15_calcw (MaxWater, maxcell),     &
                 pairlist%num_nb15_calcww(MaxWater, maxcell),     &
                 pairlist%num_nb15_water (MaxWater, maxcell),     &
                 pairlist%nb15_cell      (maxcell),               &
                 pairlist%nb15_cells     (maxcell),               &
                 pairlist%nb15_cellw     (maxcell),               &
                 pairlist%nb15_cellww    (maxcell),               &
                 pairlist%nb15_list      (MaxAtom, maxcell),      &
                 pairlist%nb15_lists     (MaxAtom, maxcell),      &
                 pairlist%nb15_listw     (MaxWater, maxcell),     &
                 pairlist%nb15_listww    (MaxWater, maxcell),     &
                 pairlist%nb15_calc_list1(MaxNb15, var_size),     &
                 pairlist%nb15_calc_list (MaxNb15, maxcell),      &
                 pairlist%nb15_calc_water(MaxNb15Water, maxcell), &
                 stat=alloc_stat)

      pairlist%num_nb15_calc1 (1:MaxAtom,  1:var_size)     = 0
      pairlist%num_nb15_calc  (1:MaxAtom,  1:maxcell)      = 0
      pairlist%num_nb15_calcs (1:MaxAtom,  1:maxcell)      = 0
      pairlist%num_nb15_calcw (1:MaxWater, 1:maxcell)      = 0
      pairlist%num_nb15_calcww(1:MaxWater, 1:maxcell)      = 0
      pairlist%num_nb15_water (1:MaxWater, 1:maxcell)      = 0

    case (PairListNoTable_CG)

#ifdef KCOMP
      if (allocated(pairlist%num_nb15_calc1)) then
        if (size(pairlist%num_nb15_calc1) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1,  &
                     pairlist%num_nb15_calc,   &
                     pairlist%nb15_calc_list1, &
                     pairlist%nb15_calc_list,  &
                     pairlist%nb15_cell,       &
                     pairlist%nb15_list,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1)) &
        allocate(pairlist%num_nb15_calc1 (MaxAtom, var_size), &
                 pairlist%num_nb15_calc  (MaxAtom, maxcell),  &
                 pairlist%nb15_calc_list1(MaxNb15_CG, var_size), &
                 pairlist%nb15_calc_list (MaxNb15_CG, maxcell),  &
                 pairlist%nb15_cell      (maxcell),           &
                 pairlist%nb15_list      (MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%num_nb15_calc1 (1:MaxAtom, 1:var_size)  = 0
      pairlist%num_nb15_calc  (1:MaxAtom, 1:maxcell)   = 0
      MaxNb15_chk = MaxNB15_CG 
#else
      call error_msg('Alloc_Pairlist> Only K computer is available')
#endif

    end select

    if (alloc_stat /=0)   call error_msg_alloc
    if (dealloc_stat /=0) call error_msg_dealloc

    return

  end subroutine alloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_pairlist
  !> @brief        deallocate pairlist information
  !! @authors      JJ
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

    ! deallocate selected variables
    !
    select case (variable)

    case (PairListNoTable)

      if (allocated(pairlist%num_nb15_calc1)) then
        deallocate(pairlist%num_nb15_calc1,  &
                   pairlist%num_nb15_calc,   &
                   pairlist%nb15_calc_list1, &
                   pairlist%nb15_calc_list,  &
                   pairlist%nb15_cell,       &
                   pairlist%nb15_list,       &
                   stat = dealloc_stat)
      end if

    case (PairListTable)

      if (allocated(pairlist%num_nb15_calc1)) then
        deallocate(pairlist%num_nb15_calc1,  &
                   pairlist%num_nb15_calc,   &
                   pairlist%num_nb15_calcs,  &
                   pairlist%num_nb15_calcw,  &
                   pairlist%num_nb15_calcww, &
                   pairlist%num_nb15_water,  &
                   pairlist%nb15_cell,       &
                   pairlist%nb15_cells,      &
                   pairlist%nb15_cellw,      &
                   pairlist%nb15_cellww,     &
                   pairlist%nb15_list,       &
                   pairlist%nb15_lists,      &
                   pairlist%nb15_listw,      &
                   pairlist%nb15_listww,     &
                   pairlist%nb15_calc_list1, &
                   pairlist%nb15_calc_list,  &
                   pairlist%nb15_calc_water, &
                   stat = dealloc_stat)
      end if

    end select

    if (dealloc_stat /=0) call error_msg_dealloc

    return

  end subroutine dealloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist_all
  !> @brief        deallocate all pairlist information
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist_all(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    call dealloc_pairlist(pairlist, PairListNoTable)
    call dealloc_pairlist(pairlist, PairListTable)

    return

  end subroutine dealloc_pairlist_all

end module sp_pairlist_str_mod
