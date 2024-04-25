!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   sp_domain_str
!> @brief   structure of domain
!! @authors Jaewoon Jung (JJ) 
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_domain_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_domain_water

    integer                       :: atom_cls_no(4)
    real(wp)                      :: charge(4)
    real(wp)                      :: mass(4)
    integer,          allocatable :: move(:)
    integer,          allocatable :: stay(:)
    integer,          allocatable :: move_integer(:,:,:)
    integer,          allocatable :: stay_integer(:,:,:)
    real(dp),         allocatable :: move_real(:,:,:)
    real(dp),         allocatable :: stay_real(:,:,:)

  end type s_domain_water

  type, public :: s_domain

    type(s_domain_water)          :: water

    integer                       :: num_atom_all
    integer                       :: num_deg_freedom
    integer                       :: num_cell_local
    integer                       :: num_cell_boundary
    integer                       :: max_num_atom

    integer                       :: cell_start(3)
    integer                       :: cell_end(3)
    integer                       :: cell_length(3)
    integer                       :: iproc_lower(3)
    integer                       :: iproc_upper(3)
    integer                       :: neighbor(-1:1,-1:1,-1:1)

    ! system size
    real(wp)                      :: system_size(3)
    real(wp)                      :: cell_size(3)

    !  DomainCellGlobal
    integer,          allocatable :: cell_g2l(:)
    integer,          allocatable :: cell_g2b(:)
    integer,          allocatable :: cell_gxyz2l(:,:,:)
    !  DomainCellLocal
    integer,          allocatable :: cell_l2g(:)
    !  DomainCellLocalBoundary
    integer,          allocatable :: cell_l2gx(:)
    integer,          allocatable :: cell_l2gy(:)
    integer,          allocatable :: cell_l2gz(:)
    integer,          allocatable :: cell_l2gx_orig(:)
    integer,          allocatable :: cell_l2gy_orig(:)
    integer,          allocatable :: cell_l2gz_orig(:)
    !  DomainCellBoundary
    integer,          allocatable :: cell_b2g(:)
    !  DomainCellPair
    integer,          allocatable :: cell_pair(:,:)
    integer,          allocatable :: virial_check(:,:)
#ifndef PGICUDA
    real(wp),         allocatable :: cell_move(:,:,:)
#else
    real(wp), allocatable, pinned :: cell_move(:,:,:)
#endif
    !  DomainCellPairList
#ifndef PGICUDA
    integer,          allocatable :: cell_pairlist1(:,:)
    integer,          allocatable :: cell_pairlist2(:,:)
    integer,          allocatable :: univ_cell_pairlist1(:,:)  ! for GPU
    integer,          allocatable :: univ_cell_pairlist2(:,:)  ! ...
#else
    integer,  allocatable, pinned :: cell_pairlist1(:,:)
    integer,  allocatable, pinned :: cell_pairlist2(:,:)
    integer,  allocatable, pinned :: univ_cell_pairlist1(:,:)  ! for GPU
    integer,  allocatable, pinned :: univ_cell_pairlist2(:,:)  ! ...
#endif
    !  DomainNeighbourCell
    integer,          allocatable :: near_neighbor_cell_count(:)
    integer,          allocatable :: neighbor_cell_count(:)
    integer,          allocatable :: neighbor_cells(:,:)
    !  DomainDynvar
    integer,          allocatable :: id_l2g(:,:)
#ifndef PGICUDA
    integer,          allocatable :: num_atom(:)
#else
    integer,  allocatable, pinned :: num_atom(:)
#endif
    integer,          allocatable :: num_atom_t0(:)
    integer,          allocatable :: num_solute(:)
    integer,          allocatable :: num_water(:)
    integer,          allocatable :: solute_list(:,:)
    integer,          allocatable :: water_list(:,:,:)
#ifndef PGICUDA
    integer,          allocatable :: atom_cls_no(:,:)
    real(wp),         allocatable :: trans_vec(:,:,:)
    real(wp),         allocatable :: translated(:,:,:)
    real(wp),         allocatable :: charge(:,:)
    real(wp),         allocatable :: force_omp(:,:,:,:)
    real(wp),         allocatable :: force_pbc(:,:,:,:)
    real(dp),         allocatable :: coord(:,:,:)
#else
    integer,  allocatable, pinned :: atom_cls_no(:,:)
    real(wp), allocatable, pinned :: trans_vec(:,:,:)    ! often used as trans1(:,:,:)
    real(wp), allocatable, pinned :: translated(:,:,:)   ! often used as trans2(:,:,:)
    real(wp), allocatable, pinned :: charge(:,:)
    real(wp), allocatable, pinned :: force_omp(:,:,:,:)
    real(wp), allocatable, pinned :: force_pbc(:,:,:,:)
    real(dp), allocatable, pinned :: coord(:,:,:)
#endif
    real(dp),         allocatable :: coord_ref(:,:,:)
    real(dp),         allocatable :: coord_old(:,:,:)
    real(dp),         allocatable :: velocity(:,:,:)
    real(dp),         allocatable :: velocity_ref(:,:,:)
    real(dp),         allocatable :: velocity_full(:,:,:)
    real(dp),         allocatable :: force(:,:,:)
    real(dp),         allocatable :: force_short(:,:,:)
    real(dp),         allocatable :: force_long(:,:,:)
    real(dp),         allocatable :: mass(:,:)
    real(dp),         allocatable :: random(:)
    real(dp),         allocatable :: virial_cellpair(:,:)
    !  DomainGlobal
    integer,          allocatable :: id_g2l(:,:)
    !  DomainPtlMove
    real(dp),         allocatable :: buf_real(:,:,:)
    integer,          allocatable :: buf_integer(:,:,:)
    integer,          allocatable :: ptl_add(:)
    integer,          allocatable :: ptl_exit(:)
    integer,          allocatable :: ptl_exit_index(:,:)

    ! FEP
    logical                       :: fep_use = .false. ! FEP flag
    integer                       :: num_atom_single_all
    real(wp)                      :: lambljA   ! lambda for LJ of part A
    real(wp)                      :: lambljB   ! lambda for LJ of part B
    real(wp)                      :: lambelA   ! lambda for elec of part A
    real(wp)                      :: lambelB   ! lambda for elec of part B
    real(wp)                      :: lambbondA ! lambda for bond of single A
    real(wp)                      :: lambbondB ! lambda for bond of single B
    real(wp)                      :: lambrest  ! lambda for restraint
    ! number of preserved atoms
    integer,          allocatable :: num_atom_preserve(:)
    ! number of preserved atoms and appearing atoms
    integer,          allocatable :: num_atom_appear_gr(:)
    ! number of preserved atoms and vanishing atoms
    integer,          allocatable :: num_atom_vanish_gr(:)
    ! list of preserved atoms for PME calculations
    integer,          allocatable :: pmelist_preserve(:,:)
    ! list of preserved atoms and appearing atoms for PME calculations
    integer,          allocatable :: pmelist_appear_gr(:,:)
    ! list of preserved atoms and vanishing atoms for PME calculations
    integer,          allocatable :: pmelist_vanish_gr(:,:)
    ! number of single region. These are used to match coordinates
    ! and velocities of single A and single B.
    integer,          allocatable :: num_atom_singleA(:)
    integer,          allocatable :: num_atom_singleB(:)
    ! id_single?(i, j, 1) is the atom index in cell j in single region.
    ! id_single?(i, j, 2) is the global atom index in single region.
    ! id_single?(i, j, 3) is the sorted index in cell j in single region.
    ! i is the index from 1 to num_atom_singleA(j).
    ! j is the cell index.
    integer,          allocatable :: id_singleA(:,:,:)
    integer,          allocatable :: id_singleB(:,:,:)
    ! temporary for forces from restraints 
    real(wp),         allocatable :: f_fep_omp(:,:,:,:)
    ! fepgrp represents atom group in FEP calculation.
    ! fepgrp has a value of 1, 2, 3, 4, or 5, which represents
    ! singleA, singleB, dualA, dualB, and preserved, respectively.
#ifndef PGICUDA
    integer,          allocatable :: fepgrp(:,:)
#else
    integer,   allocatable,pinned :: fepgrp(:,:)
#endif

  end type s_domain

  ! parameters for allocatable variables
  integer,      public, parameter :: DomainCellGlobal       = 1
  integer,      public, parameter :: DomainCellLocal        = 2
  integer,      public, parameter :: DomainCellLocBou       = 3
  integer,      public, parameter :: DomainCellBoundary     = 4
  integer,      public, parameter :: DomainCellPair         = 5
  integer,      public, parameter :: DomainCellPairList     = 6
  integer,      public, parameter :: DomainNeighbourCell    = 7
  integer,      public, parameter :: DomainDynvar           = 8
  integer,      public, parameter :: DomainDynvar_Atom      = 9
  integer,      public, parameter :: DomainGlobal           = 10
  integer,      public, parameter :: DomainPtlMove          = 11
  integer,      public, parameter :: DomainWaterMove        = 12
  integer,      public, parameter :: DomainUnivCellPairList = 13
  integer,      public, parameter :: DomainDynvar_FEP       = 14
  integer,      public, parameter :: DomainDynvar_Atom_FEP  = 15

  ! variables for maximum numbers in one cell
  integer,      public            :: MaxAtom                = 150
  integer,      public            :: MaxWater               = 50
  integer,      public            :: MaxMove                = 30
  integer,      public            :: MaxWaterMove           = 20

  ! variables for maximum cells
  integer,      public            :: maxcell, maxcell_near
  integer,      public            :: univ_maxcell, univ_maxcell1  ! for GPU
  integer,      public            :: univ_ncell_near  ! for GPU
  integer,      public            :: univ_natom_max   ! for GPU
  integer,      public            :: univ_gpu_start   ! for GPU

  ! variables : cpu time of each processor
  real(dp),     public            :: calc_time
  real(dp),     public            :: calc_time_prev

  ! subroutines
  public  :: init_domain
  public  :: alloc_domain
  public  :: dealloc_domain
  public  :: dealloc_domain_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_domain
  !> @brief        initialize domain information
  !! @authors      JJ
  !! @param[out]   domain  : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_domain(domain)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain


    domain%water%atom_cls_no(1:4)   = 0
    domain%water%charge(1:4)        = 0.0_wp
    domain%water%mass(1:4)          = 0.0_wp
    domain%num_deg_freedom          = 0
    domain%num_cell_local           = 0
    domain%num_cell_boundary        = 0
    domain%max_num_atom             = 0
    domain%num_atom_all             = 0
    domain%cell_start(1:3)          = 0
    domain%cell_end(1:3)            = 0
    domain%cell_length(1:3)         = 0
    domain%iproc_lower(1:3)         = 0
    domain%iproc_upper(1:3)         = 0
    domain%neighbor(-1:1,-1:1,-1:1) = 0

    return

  end subroutine init_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_domain
  !> @brief        allocate domain information
  !! @authors      JJ    
  !! @param[inout] domain    : domain information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable
  !! @param[in]    var_size2 : 3rd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_domain(domain, variable, var_size, var_size1, var_size2)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (DomainCellGlobal)

      if (allocated(domain%cell_g2l)) then
        if (size(domain%cell_g2l(:)) /= var_size*var_size1*var_size2) &
          deallocate(domain%cell_g2l,    &
                     domain%cell_g2b,    &
                     domain%cell_gxyz2l, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_g2l)) &
        allocate(domain%cell_g2l   (var_size*var_size1*var_size2),             &
                 domain%cell_g2b   (var_size*var_size1*var_size2),             &
                 domain%cell_gxyz2l(0:var_size+1,0:var_size1+1,0:var_size2+1), &
                 stat = alloc_stat)

      domain%cell_g2l   (1:var_size*var_size1*var_size2) = 0
      domain%cell_g2b   (1:var_size*var_size1*var_size2) = 0
      domain%cell_gxyz2l(0:var_size+1, 0:var_size1+1, 0:var_size2+1) = 0

    case (DomainCellLocal)

      if (allocated(domain%cell_l2g)) then
        if (size(domain%cell_l2g(:)) /= var_size) &
          deallocate(domain%cell_l2g, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_l2g)) &
        allocate(domain%cell_l2g(var_size), stat = alloc_stat)

      domain%cell_l2g(1:var_size) = 0

    case (DomainCellLocBou)

      if (allocated(domain%cell_l2gx)) then
        if (size(domain%cell_l2gx(:)) /= var_size) &
          deallocate(domain%cell_l2gx,      &
                     domain%cell_l2gy,      &
                     domain%cell_l2gz,      &
                     domain%cell_l2gx_orig, &
                     domain%cell_l2gy_orig, &
                     domain%cell_l2gz_orig, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_l2gx))      &
        allocate(domain%cell_l2gx     (var_size), &
                 domain%cell_l2gy     (var_size), &
                 domain%cell_l2gz     (var_size), &
                 domain%cell_l2gx_orig(var_size), &
                 domain%cell_l2gy_orig(var_size), &
                 domain%cell_l2gz_orig(var_size), &
                 stat = alloc_stat)

      domain%cell_l2gx     (1:var_size) = 0
      domain%cell_l2gy     (1:var_size) = 0
      domain%cell_l2gz     (1:var_size) = 0
      domain%cell_l2gx_orig(1:var_size) = 0
      domain%cell_l2gy_orig(1:var_size) = 0
      domain%cell_l2gz_orig(1:var_size) = 0

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        if (size(domain%cell_b2g(:)) /= var_size) &
          deallocate(domain%cell_b2g, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_b2g)) &
        allocate(domain%cell_b2g(var_size), stat = alloc_stat)

      domain%cell_b2g(1:var_size) = 0

    case (DomainCellPair)

      if (allocated(domain%cell_pair)) then
        if (size(domain%cell_pair(1,:)) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%cell_move)
#endif
          deallocate(domain%cell_pair,    &
                     domain%virial_check, &
                     domain%cell_move,    &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%cell_pair)) then
        allocate(domain%cell_pair   (   var_size, var_size), &
                 domain%virial_check(   var_size, var_size), &
                 domain%cell_move   (3, var_size, var_size), &
                 stat = alloc_stat)
#ifdef USE_GPU
        if (wp == sp) then
          call set_pinned_memory(domain%cell_move, 3*var_size*var_size*4)
        else
          call set_pinned_memory(domain%cell_move, 3*var_size*var_size*8)
        end if
#endif
      end if

      domain%cell_pair   (     1:var_size, 1:var_size) = 0
      domain%virial_check(     1:var_size, 1:var_size) = 1
      domain%cell_move   (1:3, 1:var_size, 1:var_size) = 0.0_wp

    case (DomainCellPairList)

      if (allocated(domain%cell_pairlist1)) then
        if (size(domain%cell_pairlist1(1,:)) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%cell_pairlist1)
          call unset_pinned_memory(domain%cell_pairlist2)
#endif
          deallocate(domain%cell_pairlist1,  &
                     domain%cell_pairlist2,  &
                     domain%virial_cellpair, &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%cell_pairlist1)) then
        allocate(domain%cell_pairlist1(2, var_size),          &
                 domain%cell_pairlist2(var_size1, var_size1), &
                 domain%virial_cellpair(3, var_size),         &
                 stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%cell_pairlist1, 2*var_size*4)
        call set_pinned_memory(domain%cell_pairlist2, var_size1*var_size1*4)
#endif
      end if
      domain%cell_pairlist1(1:2, 1:var_size)          = 0
      domain%cell_pairlist2(1:var_size1, 1:var_size1) = 0

    case (DomainUnivCellPairList)

      ! for GPU
      if (allocated(domain%univ_cell_pairlist1)) then
        if (size(domain%univ_cell_pairlist1(1,:)) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%univ_cell_pairlist1)
          call unset_pinned_memory(domain%univ_cell_pairlist2)
#endif
          deallocate(domain%univ_cell_pairlist1, &
                     domain%univ_cell_pairlist2, &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%univ_cell_pairlist1)) then
        allocate(domain%univ_cell_pairlist1(2, var_size),          &
                 domain%univ_cell_pairlist2(var_size1, var_size1), &
                 stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%univ_cell_pairlist1, 2*var_size*4)
        call set_pinned_memory(domain%univ_cell_pairlist2, var_size1*var_size1*4)
#endif
      end if
      domain%univ_cell_pairlist1(1:2, 1:var_size)          = 0
      domain%univ_cell_pairlist2(1:var_size1, 1:var_size1) = 0

    case (DomainNeighbourCell)

      if (allocated(domain%neighbor_cell_count)) then
        if (size(domain%neighbor_cell_count(:)) /= var_size) &
          deallocate(domain%neighbor_cell_count,      &
                     domain%near_neighbor_cell_count, &
                     domain%neighbor_cells,           &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%neighbor_cell_count)) &
        allocate(domain%neighbor_cell_count(var_size),      &
                 domain%near_neighbor_cell_count(var_size), &
                 domain%neighbor_cells(125, var_size),      &
                 stat = alloc_stat)

      domain%neighbor_cell_count  (1:var_size)      = 0
      domain%near_neighbor_cell_count  (1:var_size) = 0
      domain%neighbor_cells(1:125, 1:var_size)      = 0

    case (DomainDynvar)

      if (allocated(domain%num_atom)) then
        if (size(domain%num_atom) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%num_atom)
#endif
          deallocate(domain%num_atom,          &
                     domain%num_atom_t0,       &
                     domain%num_solute,        &
                     domain%num_water,         &
                     domain%random,            &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%num_atom)) then
          allocate(domain%num_atom(var_size),    &
                   domain%num_atom_t0(var_size), &
                   domain%num_solute(var_size),  &
                   domain%num_water(var_size),   &
                   domain%random(var_size),      &
                   stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%num_atom, var_size*4)
#endif
      end if
      domain%num_atom   (1:var_size) = 0
      domain%num_atom_t0(1:var_size) = 0
      domain%num_solute (1:var_size) = 0
      domain%num_water  (1:var_size) = 0
      domain%random     (1:var_size) = 0.0_wp

    case (DomainDynvar_FEP)

      if (allocated(domain%num_atom)) then
        if (size(domain%num_atom) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%num_atom)
#endif
          deallocate(domain%num_atom,          &
                     domain%num_atom_preserve,  &
                     domain%num_atom_singleA,   &
                     domain%num_atom_singleB,   &
                     domain%num_atom_appear_gr, &
                     domain%num_atom_vanish_gr, &
                     domain%num_atom_t0,       &
                     domain%num_solute,        &
                     domain%num_water,         &
                     domain%random,            &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%num_atom)) then
          allocate(domain%num_atom(var_size),    &
                   domain%num_atom_preserve(var_size),    &
                   domain%num_atom_singleA(var_size),     &
                   domain%num_atom_singleB(var_size),     &
                   domain%num_atom_appear_gr(var_size),   &
                   domain%num_atom_vanish_gr(var_size),   &
                   domain%num_atom_t0(var_size), &
                   domain%num_solute(var_size),  &
                   domain%num_water(var_size),   &
                   domain%random(var_size),      &
                   stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%num_atom, var_size*4)
#endif
      end if
      domain%num_atom   (1:var_size) = 0
      domain%num_atom_preserve   (1:var_size) = 0
      domain%num_atom_singleA    (1:var_size) = 0
      domain%num_atom_singleB    (1:var_size) = 0
      domain%num_atom_appear_gr  (1:var_size) = 0
      domain%num_atom_vanish_gr  (1:var_size) = 0
      domain%num_atom_t0(1:var_size) = 0
      domain%num_solute (1:var_size) = 0
      domain%num_water  (1:var_size) = 0
      domain%random     (1:var_size) = 0.0_wp

    case (DomainDynvar_Atom)

      if (allocated(domain%id_l2g)) then
        if (size(domain%id_l2g(MaxAtom,:)) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%atom_cls_no)
          call unset_pinned_memory(domain%coord)
          call unset_pinned_memory(domain%trans_vec)
          call unset_pinned_memory(domain%translated)
          call unset_pinned_memory(domain%force_omp)
          call unset_pinned_memory(domain%force_pbc)
          call unset_pinned_memory(domain%charge)
#endif
          deallocate(domain%id_l2g,        &
                     domain%solute_list,   &
                     domain%water_list,    &
                     domain%atom_cls_no,   &
                     domain%coord,         &
                     domain%coord_ref,     &
                     domain%coord_old,     &
                     domain%velocity,      &
                     domain%velocity_ref,  &
                     domain%velocity_full, &
                     domain%trans_vec,     &
                     domain%translated,    &
                     domain%force,         &
                     domain%force_short,   &
                     domain%force_long,    &
                     domain%force_omp,     &
                     domain%force_pbc,     &
                     domain%charge,        &
                     domain%mass,          &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%id_l2g)) then
        allocate(domain%id_l2g       (MaxAtom, var_size),     &
                 domain%solute_list  (MaxAtom, var_size),     &
                 domain%water_list   (var_size1, MaxWater, var_size), &
                 domain%atom_cls_no  (MaxAtom, var_size),     &
                 domain%coord        (3, MaxAtom, var_size),  &
                 domain%coord_ref    (3, MaxAtom, var_size),  &
                 domain%coord_old    (3, MaxAtom, var_size),  &
                 domain%velocity     (3, MaxAtom, var_size),  &
                 domain%velocity_ref (3, MaxAtom, var_size),  &
                 domain%velocity_full(3, MaxAtom, var_size),  &
                 domain%trans_vec    (3, MaxAtom, var_size),  &
                 domain%translated   (3, MaxAtom, var_size),  &
                 domain%force        (3, MaxAtom, var_size),  &
                 domain%force_short  (3, MaxAtom, var_size),  &
                 domain%force_long   (3, MaxAtom, var_size),  &
                 domain%force_omp    (3, MaxAtom, var_size, nthread), &
                 domain%force_pbc    (3, MaxAtom, var_size, nthread), &
                 domain%charge       (MaxAtom, var_size),     &
                 domain%mass         (MaxAtom, var_size),     &
                 stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%atom_cls_no, MaxAtom*var_size*4)
        call set_pinned_memory(domain%coord, 3*MaxAtom*var_size*8)
        if (wp == sp) then
          call set_pinned_memory(domain%trans_vec, 3*MaxAtom*var_size*4)
          call set_pinned_memory(domain%translated,3*MaxAtom*var_size*4)
          call set_pinned_memory(domain%force_omp, 3*MaxAtom*var_size*nthread*4)
          call set_pinned_memory(domain%force_pbc, 3*MaxAtom*var_size*nthread*4)
          call set_pinned_memory(domain%charge,      MaxAtom*var_size*4)
        else
          call set_pinned_memory(domain%trans_vec, 3*MaxAtom*var_size*8)
          call set_pinned_memory(domain%translated,3*MaxAtom*var_size*8)
          call set_pinned_memory(domain%force_omp, 3*MaxAtom*var_size*nthread*8)
          call set_pinned_memory(domain%force_pbc, 3*MaxAtom*var_size*nthread*8)
          call set_pinned_memory(domain%charge,      MaxAtom*var_size*8)
        end if
#endif
      end if

      domain%id_l2g       (1:MaxAtom, 1:var_size)        = 0
      domain%solute_list  (1:MaxAtom, 1:var_size)        = 0
      domain%water_list   (1:var_size1, 1:MaxWater, 1:var_size ) = 0
      domain%atom_cls_no  (1:MaxAtom, 1:var_size)        = 0
      domain%coord        (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%coord_ref    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%coord_old    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity     (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity_ref (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity_full(1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%trans_vec    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%translated   (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force        (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%force_short  (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force_long   (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force_omp    (1:3, 1:MaxAtom, 1:var_size, 1:nthread) = 0.0_wp
      domain%force_pbc    (1:3, 1:MaxAtom, 1:var_size, 1:nthread) = 0.0_wp
      domain%charge       (1:MaxAtom, 1:var_size)        = 0.0_wp
      domain%mass         (1:MaxAtom, 1:var_size)        = 0.0_wp

    case (DomainDynvar_Atom_FEP)

      if (allocated(domain%id_l2g)) then
        if (size(domain%id_l2g(MaxAtom,:)) /= var_size) then
#ifdef USE_GPU
          call unset_pinned_memory(domain%atom_cls_no)
          call unset_pinned_memory(domain%coord)
          call unset_pinned_memory(domain%trans_vec)
          call unset_pinned_memory(domain%translated)
          call unset_pinned_memory(domain%force_omp)
          call unset_pinned_memory(domain%force_pbc)
          call unset_pinned_memory(domain%charge)
          call unset_pinned_memory(domain%fepgrp)
#endif
          deallocate(domain%id_l2g,        &
                     domain%solute_list,   &
                     domain%water_list,    &
                     domain%atom_cls_no,   &
                     domain%coord,         &
                     domain%coord_ref,     &
                     domain%coord_old,     &
                     domain%velocity,      &
                     domain%velocity_ref,  &
                     domain%velocity_full, &
                     domain%trans_vec,     &
                     domain%translated,    &
                     domain%force,         &
                     domain%force_short,   &
                     domain%force_long,    &
                     domain%force_omp,     &
                     domain%force_pbc,     &
                     domain%charge,        &
                     domain%mass,          &
                     domain%pmelist_preserve,    &
                     domain%pmelist_appear_gr,   &
                     domain%pmelist_vanish_gr,   &
                     domain%id_singleA,     &
                     domain%id_singleB,     &
                     domain%fepgrp,             &
                     domain%f_fep_omp,      &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%id_l2g)) then
        allocate(domain%id_l2g       (MaxAtom, var_size),     &
                 domain%solute_list  (MaxAtom, var_size),     &
                 domain%water_list   (var_size1, MaxWater, var_size), &
                 domain%atom_cls_no  (MaxAtom, var_size),     &
                 domain%coord        (3, MaxAtom, var_size),  &
                 domain%coord_ref    (3, MaxAtom, var_size),  &
                 domain%coord_old    (3, MaxAtom, var_size),  &
                 domain%velocity     (3, MaxAtom, var_size),  &
                 domain%velocity_ref (3, MaxAtom, var_size),  &
                 domain%velocity_full(3, MaxAtom, var_size),  &
                 domain%trans_vec    (3, MaxAtom, var_size),  &
                 domain%translated   (3, MaxAtom, var_size),  &
                 domain%force        (3, MaxAtom, var_size),  &
                 domain%force_short  (3, MaxAtom, var_size),  &
                 domain%force_long   (3, MaxAtom, var_size),  &
                 domain%force_omp    (3, MaxAtom, var_size, nthread), &
                 domain%force_pbc    (3, MaxAtom, var_size, nthread), &
                 domain%charge       (MaxAtom, var_size),     &
                 domain%mass         (MaxAtom, var_size),     &
                 domain%pmelist_preserve  (MaxAtom, var_size),  &
                 domain%pmelist_appear_gr (MaxAtom, var_size),  &
                 domain%pmelist_vanish_gr (MaxAtom, var_size),  &
                 domain%id_singleA     (MaxAtom, var_size, 3),  &
                 domain%id_singleB     (MaxAtom, var_size, 3),  &
                 domain%fepgrp            (MaxAtom, var_size),  &
                 domain%f_fep_omp(3, MaxAtom, var_size, nthread), &
                 stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(domain%atom_cls_no, MaxAtom*var_size*4)
        call set_pinned_memory(domain%coord, 3*MaxAtom*var_size*8)
        if (wp == sp) then
          call set_pinned_memory(domain%trans_vec, 3*MaxAtom*var_size*4)
          call set_pinned_memory(domain%translated,3*MaxAtom*var_size*4)
          call set_pinned_memory(domain%force_omp, 3*MaxAtom*var_size*nthread*4)
          call set_pinned_memory(domain%force_pbc, 3*MaxAtom*var_size*nthread*4)
          call set_pinned_memory(domain%charge,      MaxAtom*var_size*4)
        else
          call set_pinned_memory(domain%trans_vec, 3*MaxAtom*var_size*8)
          call set_pinned_memory(domain%translated,3*MaxAtom*var_size*8)
          call set_pinned_memory(domain%force_omp, 3*MaxAtom*var_size*nthread*8)
          call set_pinned_memory(domain%force_pbc, 3*MaxAtom*var_size*nthread*8)
          call set_pinned_memory(domain%charge,      MaxAtom*var_size*8)
        end if
        call set_pinned_memory(domain%fepgrp, MaxAtom*var_size*4)
#endif
      end if

      domain%id_l2g       (1:MaxAtom, 1:var_size)        = 0
      domain%solute_list  (1:MaxAtom, 1:var_size)        = 0
      domain%water_list   (1:var_size1, 1:MaxWater, 1:var_size ) = 0
      domain%atom_cls_no  (1:MaxAtom, 1:var_size)        = 0
      domain%coord        (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%coord_ref    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%coord_old    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity     (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity_ref (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%velocity_full(1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%trans_vec    (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%translated   (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force        (1:3, 1:MaxAtom, 1:var_size)   = 0.0_dp
      domain%force_short  (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force_long   (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp
      domain%force_omp    (1:3, 1:MaxAtom, 1:var_size, 1:nthread) = 0.0_wp
      domain%force_pbc    (1:3, 1:MaxAtom, 1:var_size, 1:nthread) = 0.0_wp
      domain%charge       (1:MaxAtom, 1:var_size)        = 0.0_wp
      domain%mass         (1:MaxAtom, 1:var_size)        = 0.0_wp
      domain%pmelist_preserve  (1:MaxAtom, 1:var_size)   = 0
      domain%pmelist_appear_gr (1:MaxAtom, 1:var_size)   = 0
      domain%pmelist_vanish_gr (1:MaxAtom, 1:var_size)   = 0
      domain%id_singleA   (1:MaxAtom, 1:var_size, 1:3)   = 0
      domain%id_singleB   (1:MaxAtom, 1:var_size, 1:3)   = 0
      domain%fepgrp            (1:MaxAtom, 1:var_size)   = 5
      domain%f_fep_omp(1:3, 1:MaxAtom, 1:var_size, 1:nthread) = 0.0_wp

    case (DomainGlobal)

      if (allocated(domain%id_g2l)) then
        if (size(domain%id_g2l(2,:)) /= var_size) &
          deallocate(domain%id_g2l, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_g2l)) &
        allocate(domain%id_g2l(2,var_size), stat = alloc_stat)

      domain%id_g2l(1:2,1:var_size) = 0

    case (DomainPtlMove) 

      if (allocated(domain%buf_real)) then
        if (size(domain%buf_real(1,1,:)) /= var_size) &
          deallocate(domain%buf_real,       &
                     domain%buf_integer,    &
                     domain%ptl_add,        &
                     domain%ptl_exit,       &
                     domain%ptl_exit_index, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%buf_real)) &
        allocate(domain%buf_real      (8, MaxMove, var_size), &
                 domain%buf_integer   (9, MaxMove, var_size), &
                 domain%ptl_add       (var_size),             &
                 domain%ptl_exit      (var_size1),            &
                 domain%ptl_exit_index(MaxMove, var_size1),   &
                 stat = alloc_stat)

      domain%buf_real      (1:8, 1:MaxMove, 1:var_size) = 0.0_dp
      domain%buf_integer   (1:9, 1:MaxMove, 1:var_size) = 0
      domain%ptl_add       (1:var_size)                 = 0
      domain%ptl_exit      (1:var_size1)                = 0
      domain%ptl_exit_index(1:MaxMove, 1:var_size1)     = 0

    case (DomainWaterMove)

      if (allocated(domain%water%move)) then
        if (size(domain%water%move(:)) /= var_size) &
          deallocate(domain%water%move,         &
                     domain%water%stay,         &
                     domain%water%move_integer, &
                     domain%water%stay_integer, &
                     domain%water%move_real,    &
                     domain%water%stay_real,    &
                     stat = dealloc_stat )
      end if

      if (.not. allocated(domain%water%move)) &
        allocate(domain%water%move        (var_size),                            &
                 domain%water%stay        (var_size1),                           &
                 domain%water%move_integer(2*var_size2, MaxWaterMove, var_size),   &
                 domain%water%stay_integer(3*var_size2, MaxWater, var_size1),      &
                 domain%water%move_real   (6*var_size2, MaxWaterMove, var_size), &
                 domain%water%stay_real   (6*var_size2, MaxWater, var_size1),    &
                 stat = alloc_stat)

      domain%water%move        (1:var_size)                       = 0
      domain%water%stay        (1:var_size1)                      = 0
      domain%water%move_integer(1:2*var_size2, 1:MaxWaterMove, 1:var_size)  = 0
      domain%water%stay_integer(1:2*var_size2, 1:MaxWater, 1:var_size1)     = 0
      domain%water%move_real   (1:6*var_size2, 1:MaxWaterMove, 1:var_size) = 0.0_dp
      domain%water%stay_real   (1:6*var_size2, 1:MaxWater, 1:var_size1)    = 0.0_dp

    case default

      call error_msg('Alloc_Domain> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain
  !> @brief        deallocate domain information
  !! @authors      JJ    
  !! @@aram[inout] domain   : domain information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_domain(domain, variable)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    select case (variable)


    case (DomainCellGlobal)

      if (allocated(domain%cell_g2l)) then
        deallocate(domain%cell_g2l,    &
                   domain%cell_g2b,    &
                   domain%cell_gxyz2l, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocal)

      if (allocated(domain%cell_l2g)) then
        deallocate(domain%cell_l2g, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocBou)

      if (allocated(domain%cell_l2gx)) then
        deallocate(domain%cell_l2gx,      &
                   domain%cell_l2gy,      &
                   domain%cell_l2gz,      &
                   domain%cell_l2gx_orig, &
                   domain%cell_l2gy_orig, &
                   domain%cell_l2gz_orig, &
                   stat = dealloc_stat)
      end if

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        deallocate(domain%cell_b2g, &
                   stat = dealloc_stat)
      end if

    case (DomainCellPair)

      if (allocated(domain%cell_pair)) then
#ifdef USE_GPU
        call unset_pinned_memory(domain%cell_move)
#endif
        deallocate(domain%cell_pair,    &
                   domain%virial_check, &
                   domain%cell_move,    &
                   stat = dealloc_stat)
      end if

    case (DomainCellPairList) 

      if (allocated(domain%cell_pairlist1)) then
        deallocate(domain%cell_pairlist1, &
                   domain%cell_pairlist2, &
                   stat = dealloc_stat)
#ifdef USE_GPU
        call unset_pinned_memory(domain%univ_cell_pairlist1)
        call unset_pinned_memory(domain%univ_cell_pairlist2)
        deallocate(domain%univ_cell_pairlist1, &  ! for GPU
                   domain%univ_cell_pairlist2, &
                   stat = dealloc_stat)
#endif
      end if

    case (DomainNeighbourCell)

      if (allocated(domain%neighbor_cell_count)) then
        deallocate(domain%neighbor_cell_count, &
                   domain%neighbor_cells,      &
                   stat = dealloc_stat)
      end if

    case (DomainDynvar)

      if (allocated(domain%coord)) then
#ifdef USE_GPU
        call unset_pinned_memory(domain%num_atom)
        call unset_pinned_memory(domain%atom_cls_no)
        call unset_pinned_memory(domain%coord)
        call unset_pinned_memory(domain%trans_vec)
        call unset_pinned_memory(domain%translated)
        call unset_pinned_memory(domain%force_omp)
        call unset_pinned_memory(domain%force_pbc)
        call unset_pinned_memory(domain%charge)
#endif
        deallocate(domain%id_l2g,        &
                   domain%num_atom,      &
                   domain%num_atom_t0,   &
                   domain%num_solute,    &
                   domain%num_water,     &
                   domain%solute_list,   &
                   domain%water_list,    &
                   domain%atom_cls_no,   &
                   domain%coord,         &
                   domain%coord_ref,     &
                   domain%coord_old,     &
                   domain%velocity,      &
                   domain%velocity_ref,  &
                   domain%velocity_full, &
                   domain%trans_vec,     &
                   domain%translated,    &
                   domain%force,         &
                   domain%force_short,   &
                   domain%force_long,    &
                   domain%force_omp,     &
                   domain%force_pbc,     &
                   domain%charge,        &
                   domain%mass,          &
                   domain%random,        &
                   stat = dealloc_stat)
      end if

    case (DomainDynvar_FEP)

      if (allocated(domain%coord)) then
#ifdef USE_GPU
        call unset_pinned_memory(domain%num_atom)
        call unset_pinned_memory(domain%atom_cls_no)
        call unset_pinned_memory(domain%coord)
        call unset_pinned_memory(domain%trans_vec)
        call unset_pinned_memory(domain%translated)
        call unset_pinned_memory(domain%force_omp)
        call unset_pinned_memory(domain%force_pbc)
        call unset_pinned_memory(domain%charge)
        call unset_pinned_memory(domain%fepgrp)
#endif
        deallocate(domain%id_l2g,                  &
                   domain%num_atom,                &
                   domain%num_atom_preserve,       &
                   domain%num_atom_singleA,        &
                   domain%num_atom_singleB,        &
                   domain%num_atom_appear_gr,      &
                   domain%num_atom_vanish_gr,      &
                   domain%pmelist_preserve,        &
                   domain%pmelist_appear_gr,       &
                   domain%pmelist_vanish_gr,       &
                   domain%id_singleA,              &
                   domain%id_singleB,              &
                   domain%num_atom_t0,             &
                   domain%num_solute,              &
                   domain%num_water,               &
                   domain%solute_list,             &
                   domain%water_list,              &
                   domain%atom_cls_no,             &
                   domain%coord,                   &
                   domain%coord_ref,               &
                   domain%coord_old,               &
                   domain%velocity,                &
                   domain%velocity_ref,            &
                   domain%velocity_full,           &
                   domain%trans_vec,               &
                   domain%translated,              &
                   domain%force,                   &
                   domain%force_short,             &
                   domain%force_long,              &
                   domain%force_omp,               &
                   domain%force_pbc,               &
                   domain%charge,                  &
                   domain%mass,                    &
                   domain%random,                  &
                   domain%fepgrp,                  &
                   domain%f_fep_omp,               &
                   stat = dealloc_stat)
      end if

    case (DomainGlobal)

      if (allocated(domain%id_g2l)) then
        deallocate(domain%id_g2l, &
                   stat = dealloc_stat)
      end if

    case (DomainPtlMove) 

      if (allocated(domain%buf_real)) then
        deallocate(domain%buf_real,       &
                   domain%buf_integer,    &
                   domain%ptl_add,        &
                   domain%ptl_exit,       &
                   domain%ptl_exit_index, &
                   stat = dealloc_stat)
      end if

    case (DomainWaterMove)

      if (allocated(domain%water%move)) then
        deallocate(domain%water%move,         &
                   domain%water%stay,         &
                   domain%water%move_integer, &
                   domain%water%stay_integer, &
                   domain%water%move_real,    &
                   domain%water%stay_real,    &
                   stat = dealloc_stat )
      end if

    case default

      call error_msg('Dealloc_Domain> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_domain


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain_all
  !> @brief        deallocate all domain information
  !! @authors      JJ
  !! @param[inout] domain : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_domain_all(domain)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain


    call dealloc_domain(domain, DomainCellGlobal)
    call dealloc_domain(domain, DomainCellLocal)
    call dealloc_domain(domain, DomainCellLocBou)
    call dealloc_domain(domain, DomainCellBoundary)
    call dealloc_domain(domain, DomainCellPair)
    call dealloc_domain(domain, DomainCellPairList)
    call dealloc_domain(domain, DomainNeighbourCell)
    call dealloc_domain(domain, DomainDynvar)
    call dealloc_domain(domain, DomainGlobal)
    call dealloc_domain(domain, DomainPtlMove)
    call dealloc_domain(domain, DomainWaterMove)

    return

  end subroutine dealloc_domain_all

end module sp_domain_str_mod
