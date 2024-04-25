!--------1---------2---------3---------4---------5---------6---------7---------8 ! 
!  Module   sp_domain_mod
!> @brief   utilities for domain decomposition                 
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_domain_mod

  use sp_constraints_mod
  use sp_boundary_mod
  use sp_energy_mod
  use sp_communicate_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use molecules_mod
  use timers_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters for setup cell capacity
  real(dp),         private, parameter :: VolumeBox8 = 512.0_dp
  real(dp),         private, parameter :: VboxRate   = 1.4_dp
  real(dp),         private, parameter :: ShrinkRate = 1.2_dp
  integer,          private, parameter :: NAtomBox8  = 120
  integer,          private, parameter :: NBondBox8  = 120
  integer,          private, parameter :: NAnglBox8  = 300
  integer,          private, parameter :: NDiheBox8  = 700
  integer,          private, parameter :: NImprBox8  = 50
  integer,          private, parameter :: NCmapBox8  = 12
  integer,          private, parameter :: NHGrpBox8  = 40
  integer,          private, parameter :: NHMovBox8  = 8
  integer,          private, parameter :: NContBox8  = 300
#ifdef USE_GPU
  integer,          private, parameter :: NAtomMax_in_CUDA = 255
#endif
#ifndef KCOMP
  real(dp),         private, parameter :: Scale_pairlist = 1.0_dp
#else
  real(dp),         private, parameter :: Scale_pairlist = 0.7_dp
#endif

  ! subroutines
  public  :: setup_domain
  public  :: setup_domain_pio
  public  :: setup_domain_interaction
  public  :: update_domain_interaction
  public  :: update_outgoing_ptl
  public  :: update_incoming_ptl
  public  :: setup_processor_rank
  public  :: setup_cell_capacity
  public  :: setup_cell_capacity_pio
  public  :: setup_cell_boundary
  public  :: update_nobc_boundary
  private :: setup_solute_and_water
  private :: setup_hbond_group
  private :: setup_atom_by_HBond
  private :: setup_atom_by_table
  private :: setup_atom
  private :: assign_neighbor_cells
  private :: assign_cell_atoms
  private :: assign_cell_cpu
  private :: assign_cell_interaction
  private :: molecule_to_domain
  private :: check_atom_coord
  private :: check_atom_coord_pio
  !FEP
  public  :: setup_domain_fep
  public  :: update_outgoing_ptl_fep
  public  :: update_incoming_ptl_fep
  public  :: setup_cell_capacity_fep
  public  :: update_pmelist_fep
  public  :: sort_single_fep
  private :: setup_hbond_group_fep
  private :: setup_atom_by_HBond_fep
  private :: setup_atom_by_table_fep
  private :: setup_atom_fep
  private :: molecule_to_domain_fep
  private :: check_atom_coord_fep
  private :: setup_pmelist_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain
  !> @brief        setup domain information
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain(ene_info, con_info, &
                          boundary, molecule, enefunc, constraints, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all


    ! initialize structure informations
    !
    call init_domain(domain)
    call init_enefunc(enefunc)
    call init_constraints(constraints)

    domain%num_atom_all       = molecule%num_atoms

    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)


    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity(boundary, domain, molecule)


    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)


    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do


    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)


    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all, 1, 1)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all,    1, 1)

    ! decide hydrogen atom from mass
    !
    call check_light_atom_name(con_info%hydrogen_mass_upper_bound, molecule)

    if (boundary%type == BoundaryTypePBC) then

      if (constraints%rigid_bond) then

        call setup_solute_and_water(molecule, enefunc,       &
                                    constraints%water_model, &
                                    constraints%tip4)

        if (constraints%tip4) then
          call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
        else
          call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
        end if

        call setup_hbond_group     (molecule, enefunc, constraints)

        call setup_atom_by_HBond   (molecule, boundary, enefunc, &
                                    constraints, domain)

      else 

        call setup_solute_and_water(molecule, enefunc,         &
                                    enefunc%table%water_model, &
                                    enefunc%table%tip4)

        if (enefunc%table%tip4) then
          call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
        else
          call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
        end if
        call setup_atom_by_table   (molecule, boundary, enefunc, domain)

      end if

    else if (boundary%type == BoundaryTypeNOBC) then

      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 0, 1)
      call setup_atom (molecule, boundary, domain)

    end if

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    call check_atom_coord(ene_info, boundary, ene_info%contact_check, domain)

    ! setup water molecule information
    !
    if (boundary%type == BoundaryTypePBC .and. &
        enefunc%table%num_water > 0) then

      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
      domain%water%charge(1)      = enefunc%table%charge_O
      domain%water%charge(2)      = enefunc%table%charge_H
      domain%water%charge(3)      = enefunc%table%charge_H
      domain%water%charge(4)      = enefunc%table%charge_D
      domain%water%mass(1)        = enefunc%table%mass_O
      domain%water%mass(2)        = enefunc%table%mass_H
      domain%water%mass(3)        = enefunc%table%mass_H
      domain%water%mass(4)        = enefunc%table%mass_D

    end if

    return

  end subroutine setup_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_pio
  !> @brief        setup domain information
  !! @authors      JJ
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_pio(con_info, &
                              boundary, enefunc, constraints, domain)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all


    ! initialize structure informations
    !
    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)


    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity_pio(boundary, domain)


    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)


    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do


    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)


    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    if (enefunc%contact_check) &
      call check_atom_coord_pio(enefunc, boundary, domain)

    ! setup water molecule information
    !
    domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
    domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
    domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
    domain%water%charge(1)      = enefunc%table%charge_O
    domain%water%charge(2)      = enefunc%table%charge_H
    domain%water%charge(3)      = enefunc%table%charge_H
    domain%water%mass(1)        = enefunc%table%mass_O
    domain%water%mass(2)        = enefunc%table%mass_H
    domain%water%mass(3)        = enefunc%table%mass_H

    return

  end subroutine setup_domain_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_interaction
  !> @brief        define the pairwise interaction between cells
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_interaction(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, j, ij, inbc
    integer                   :: ncel_local, nboundary, cell(3)

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(wp),         pointer :: cell_move(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)
    integer,          pointer :: cell_pair(:,:), cell_gxyz2l(:,:,:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx1(:), cell_l2gy1(:), cell_l2gz1(:)
    integer,          pointer :: neighbor(:,:,:), num_domain(:)
    integer,          pointer :: nsolute(:), nwater(:), natom_t0(:)
    integer,          pointer :: virial_check(:,:)
  
    real(wp),         allocatable :: natom(:)


    bsize_x     => boundary%box_size_x
    bsize_y     => boundary%box_size_y
    bsize_z     => boundary%box_size_z
    num_domain  => boundary%num_domain

    cell_start   => domain%cell_start
    cell_end     => domain%cell_end
    cell_pair    => domain%cell_pair
    cell_move    => domain%cell_move
    virial_check => domain%virial_check
    cell_gxyz2l  => domain%cell_gxyz2l
    cell_l2gx    => domain%cell_l2gx
    cell_l2gy    => domain%cell_l2gy
    cell_l2gz    => domain%cell_l2gz
    cell_l2gx1   => domain%cell_l2gx_orig
    cell_l2gy1   => domain%cell_l2gy_orig
    cell_l2gz1   => domain%cell_l2gz_orig
    nsolute      => domain%num_solute
    nwater       => domain%num_water
    natom_t0     => domain%num_atom_t0
    neighbor     => domain%neighbor

    cell(1)      = boundary%num_cells_x
    cell(2)      = boundary%num_cells_y
    cell(3)      = boundary%num_cells_z
    ncel_local   = domain%num_cell_local
    nboundary    = domain%num_cell_boundary

    allocate(natom(1:ncel_local+nboundary))

    ! check neighboring cells of each local cell
    !
    call assign_neighbor_cells(boundary, domain)

    ! number of atoms in each cell (proceesor number is also considered)
    !
    call assign_cell_atoms(natom, natom_t0, &
                           cell_l2gx, cell_l2gy, cell_l2gz, &
                           cell_start, cell_end,            &
                           neighbor, ncel_local, nboundary)

    ! assign the interaction cell for each interaction
    !
    call assign_cell_interaction(natom, &
                                 bsize_x, bsize_y, bsize_z, num_domain, cell, &
                                 cell_l2gx,  cell_l2gy,  cell_l2gz,   &
                                 cell_l2gx1, cell_l2gy1, cell_l2gz1,  &
                                 cell_gxyz2l, cell_start, cell_end,   &
                                 ncel_local, nboundary, cell_move,    &
                                 cell_pair, virial_check)

    if (boundary%type == BoundaryTypePBC) then

      ij = 0
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
          end if
        end do
      end do
      maxcell_near = ij
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1,                    &
                  domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
          end if
        end do
      end do
      maxcell = ij
      univ_maxcell = ij + ncel_local

    else if (boundary%type == BoundaryTypeNOBC) then

      ij = 0
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            if (abs(cell_move(1,j,i)) < EPS .and. &
                abs(cell_move(2,j,i)) < EPS .and. &
                abs(cell_move(3,j,i)) < EPS) then
              ij = ij + 1
            end if
          end if
        end do
      end do
      maxcell_near = ij
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1,  &
                  domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            if (abs(cell_move(1,j,i)) < EPS .and. &
                abs(cell_move(2,j,i)) < EPS .and. &
                abs(cell_move(3,j,i)) < EPS) then
              ij = ij + 1
            end if
          end if
        end do
      end do
      maxcell = ij
      univ_maxcell = ij + ncel_local

    end if
    call alloc_domain(domain, DomainCellPairList, maxcell, &
                      ncel_local+nboundary, 1)

#ifdef USE_GPU
    call alloc_domain(domain, DomainUnivCellPairList, univ_maxcell, &
                      ncel_local+nboundary, 1)
#endif

    if (boundary%type == BoundaryTypePBC) then

      ij = 0
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%cell_pairlist1(1,ij) = i
            domain%cell_pairlist1(2,ij) = j
            domain%cell_pairlist2(j,i) = ij
            domain%cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      maxcell_near = ij
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1,                       &
                  domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%cell_pairlist1(1,ij) = i
            domain%cell_pairlist1(2,ij) = j
            domain%cell_pairlist2(j,i) = ij
            domain%cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      maxcell = ij

#ifdef USE_GPU
      ij = 0

      ! cell-pairs (same cell)
      !
      do i = 1, ncel_local
         ij = ij + 1
         domain%univ_cell_pairlist1(1,ij) = i
         domain%univ_cell_pairlist1(2,ij) = i
         domain%univ_cell_pairlist2(i,i) = ij
      enddo

      ! cell-pairs (different cell)
      !
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%univ_cell_pairlist1(1,ij) = i
            domain%univ_cell_pairlist1(2,ij) = j
            domain%univ_cell_pairlist2(j,i) = ij
            domain%univ_cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1,                      &
                  domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%univ_cell_pairlist1(1,ij) = i
            domain%univ_cell_pairlist1(2,ij) = j
            domain%univ_cell_pairlist2(j,i) = ij
            domain%univ_cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      univ_maxcell = ij
      univ_ncell_near = ncel_local + maxcell_near
#endif

    else if (boundary%type == BoundaryTypeNOBC) then

      ij = 0
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            if (abs(cell_move(1,j,i)) < EPS .and. &
                abs(cell_move(2,j,i)) < EPS .and. &
                abs(cell_move(3,j,i)) < EPS) then
              ij = ij + 1
              domain%cell_pairlist1(1,ij) = i
              domain%cell_pairlist1(2,ij) = j
              domain%cell_pairlist2(j,i) = ij
              domain%cell_pairlist2(i,j) = ij
            end if
          end if
        end do
      end do
      maxcell_near = ij
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1, &
                  domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            if (abs(cell_move(1,j,i)) < EPS .and. &
                abs(cell_move(2,j,i)) < EPS .and. &
                abs(cell_move(3,j,i)) < EPS) then
              ij = ij + 1
              domain%cell_pairlist1(1,ij) = i
              domain%cell_pairlist1(2,ij) = j
              domain%cell_pairlist2(j,i) = ij
              domain%cell_pairlist2(i,j) = ij
            end if
          end if
        end do
      end do
      maxcell = ij

    end if
    
    deallocate(natom)

    return

  end subroutine setup_domain_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_domain_interaction
  !> @brief        update the pairwise interaction between cells
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_domain_interaction(boundary, enefunc, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(dp)                  :: calc_time
    integer                   :: i, j, ij, inbc
    integer                   :: ncel_local, nboundary, cell(3)

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(wp),         pointer :: cell_move(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)
    integer,          pointer :: cell_pair(:,:), cell_gxyz2l(:,:,:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx1(:), cell_l2gy1(:), cell_l2gz1(:)
    integer,          pointer :: neighbor(:,:,:), num_domain(:)
    integer,          pointer :: virial_check(:,:)

    real(wp),         allocatable :: cpu_time(:), time_proc(:)


    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    num_domain   => boundary%num_domain

    cell_start   => domain%cell_start
    cell_end     => domain%cell_end
    cell_pair    => domain%cell_pair
    cell_move    => domain%cell_move
    cell_gxyz2l  => domain%cell_gxyz2l
    cell_l2gx    => domain%cell_l2gx
    cell_l2gy    => domain%cell_l2gy
    cell_l2gz    => domain%cell_l2gz
    cell_l2gx1   => domain%cell_l2gx_orig
    cell_l2gy1   => domain%cell_l2gy_orig
    cell_l2gz1   => domain%cell_l2gz_orig
    neighbor     => domain%neighbor
    virial_check => domain%virial_check

    cell(1)      = boundary%num_cells_x
    cell(2)      = boundary%num_cells_y
    cell(3)      = boundary%num_cells_z
    ncel_local   = domain%num_cell_local
    nboundary    = domain%num_cell_boundary


    allocate(cpu_time(1:ncel_local+nboundary), time_proc(1:nproc_city))

    ! cpu time(real) of each processor
    !
    if (enefunc%pme_use) then
      calc_time = total_time(3)+total_time(4)+total_time(5)+total_time(8)       
    else      
      calc_time = total_time(3)+total_time(4)+total_time(5)+total_time(7)       
    end if

    calc_time      = calc_time - calc_time_prev
    calc_time_prev = calc_time_prev + calc_time

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(calc_time, 1, mpi_real8, time_proc, 1, &
                       mpi_real8, mpi_comm_country, ierror)
#else
    time_proc(1) = calc_time
#endif

    ! cpu time (cpu time of domain * random number of each cell)
    !
    call assign_cell_cpu(time_proc, domain%random,        &
                         cell_l2gx, cell_l2gy, cell_l2gz, &
                         cell_start, cell_end, neighbor,  &
                         ncel_local, nboundary, cpu_time)

    ! assign the interaction cell for each interaction
    !
    call assign_cell_interaction(cpu_time, &
                                 bsize_x, bsize_y, bsize_z, num_domain, cell, &
                                 cell_l2gx,  cell_l2gy,  cell_l2gz,   &
                                 cell_l2gx1, cell_l2gy1, cell_l2gz1,  &
                                 cell_gxyz2l, cell_start, cell_end,   &
                                 ncel_local, nboundary, cell_move,    &
                                 cell_pair, virial_check)

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
        end if
      end do
    end do
    maxcell = ij

    !!call alloc_domain(domain, DomainCellPairList, maxcell, &
    !!                  ncel_local+nboundary,1)
    domain%cell_pairlist1(1:2,1:maxcell) = 0
    domain%cell_pairlist2(1:ncel_local+nboundary,1:ncel_local+nboundary) = 0

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
          domain%cell_pairlist1(1,ij) = i
          domain%cell_pairlist1(2,ij) = j
          domain%cell_pairlist2(j,i) = ij
          domain%cell_pairlist2(i,j) = ij
        end if
      end do
    end do
    maxcell = ij

    deallocate(cpu_time, time_proc)

    return

  end subroutine update_domain_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_ptl
  !> @brief        update particles in each cell
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_ptl(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3)
    integer                   :: i, k, ix, icx, icy, icz, icel, ncel
    integer                   :: icel_local, icel_bd

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    real(dp),         pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(dp),         pointer :: mass(:,:)
    real(dp),         pointer :: buf_real(:,:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: atmcls(:,:), id_l2g(:,:)
    integer,          pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,          pointer :: buf_int(:,:,:)


    bsize_x  => boundary%box_size_x
    bsize_y  => boundary%box_size_y
    bsize_z  => boundary%box_size_z
    ncel_x   => boundary%num_cells_x
    ncel_y   => boundary%num_cells_y
    ncel_z   => boundary%num_cells_z
    csize_x  => boundary%cell_size_x
    csize_y  => boundary%cell_size_y
    csize_z  => boundary%cell_size_z

    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real


    ! initializaiton
    !
    ncel = domain%num_cell_local + domain%num_cell_boundary
    ptl_exit(1:domain%num_cell_local) = 0
    ptl_add (1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, domain%num_cell_local

      k = 0
      do ix = 1, domain%num_atom(i)

        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1

        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        icel_local = domain%cell_g2l(icel)
        icel_bd    = domain%cell_g2b(icel)

        if (icel_local /= i) then

          ptl_exit(i) = ptl_exit(i) + 1
          ptlindex(ptl_exit(i),i) = ix

          if (icel_local /= 0) then

            ptl_add(icel_local) = ptl_add(icel_local) + 1
            buf_real(1,ptl_add(icel_local),icel_local) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_local),icel_local) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_local),icel_local) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_local),icel_local) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_local),icel_local) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_local),icel_local) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_local),icel_local) = charge(ix,i)
            buf_real(8,ptl_add(icel_local),icel_local) = mass(ix,i)
            buf_int (1,ptl_add(icel_local),icel_local) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_local),icel_local) = id_l2g(ix,i)

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + domain%num_cell_local 
            ptl_add(icel_bd) = ptl_add(icel_bd) + 1
            buf_real(1,ptl_add(icel_bd),icel_bd) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_bd),icel_bd) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_bd),icel_bd) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_bd),icel_bd) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_bd),icel_bd) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_bd),icel_bd) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_bd),icel_bd) = charge(ix,i)
            buf_real(8,ptl_add(icel_bd),icel_bd) = mass(ix,i)
            buf_int (1,ptl_add(icel_bd),icel_bd) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_bd),icel_bd) = id_l2g(ix,i)

          end if
        end if
      end do
    end do

    return

  end subroutine update_outgoing_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_ptl
  !> @brief        update particles in each cell
  !! @authors      JJ
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_ptl(domain)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(dp),        pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: atmcls(:,:), id_l2g(:,:), id_g2l(:,:)
    integer,         pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,         pointer :: buf_int(:,:,:)


    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    id_g2l   => domain%id_g2l
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real


    ! Incoming particles
    !
    do i = 1, domain%num_cell_local

      ! When the number of coming particles is larger than that of outgoing ones
      !
      if (ptl_add(i) >= ptl_exit(i)) then

        do k = 1, ptl_exit(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        do k = ptl_exit(i)+1, ptl_add(i)
          ix = k + domain%num_atom(i) - ptl_exit(i)
          coord(1,ix,i)               = buf_real(1,k,i)
          coord(2,ix,i)               = buf_real(2,k,i)
          coord(3,ix,i)               = buf_real(3,k,i)
          velocity(1,ix,i)            = buf_real(4,k,i)
          velocity(2,ix,i)            = buf_real(5,k,i)
          velocity(3,ix,i)            = buf_real(6,k,i)
          charge(ix,i)                = buf_real(7,k,i)
          mass(ix,i)                  = buf_real(8,k,i)
          atmcls(ix,i)                = buf_int(1,k,i)
          id_l2g(ix,i)                = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ix
        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)

      ! When the number of coming particles is less than that of outgoing ones
      !
      else

        do k = 1, ptl_add(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        j  = 0
        ix = domain%num_atom(i)
        k  = ptl_add(i) + 1

        do while (j < (ptl_exit(i)-ptl_add(i)))

          insert = .true.
          do kx = k, ptl_exit(i)
            if (ix == ptlindex(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = ptlindex(k,i)
            coord(1,kx,i)          = coord(1,ix,i)
            coord(2,kx,i)          = coord(2,ix,i)
            coord(3,kx,i)          = coord(3,ix,i)
            velocity(1,kx,i)       = velocity(1,ix,i)
            velocity(2,kx,i)       = velocity(2,ix,i)
            velocity(3,kx,i)       = velocity(3,ix,i)
            charge(kx,i)           = charge(ix,i)
            mass(kx,i)             = mass(ix,i)
            atmcls(kx,i)           = atmcls(ix,i)
            id_l2g(kx,i)           = id_l2g(ix,i)
            id_g2l(1,id_l2g(kx,i)) = i
            id_g2l(2,id_l2g(kx,i)) = kx

            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)
           
      end if

      domain%num_solute(i) = domain%num_atom(i)

    end do

    return

  end subroutine update_incoming_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_rank
  !> @brief        define the processor rank in each dimension and decide the 
  !!               number of cells in each domain
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !! @@aram[out]   cell     : cells in boundary
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_rank(boundary, domain, cell)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain
    integer,                 intent(inout) :: cell(3)

    ! local variables
    integer                  :: iproc(3)
    integer                  :: i, j, k, quotient, remainder

    integer, pointer         :: ncell_local, ncell_boundary, num_domain(:)
    integer, pointer         :: cell_start(:), cell_end(:), cell_length(:)
    integer, pointer         :: iproc_lower(:), iproc_upper(:), neighbor(:,:,:)


    num_domain     => boundary%num_domain

    ncell_local    => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end  
    cell_length    => domain%cell_length
    iproc_lower    => domain%iproc_lower
    iproc_upper    => domain%iproc_upper
    neighbor       => domain%neighbor

    ! Assign the rank of each dimension from my_rank
    !
    iproc(1) = mod(my_city_rank, num_domain(1))
    iproc(2) = mod(my_city_rank/num_domain(1), num_domain(2))
    iproc(3) = my_city_rank/(num_domain(1)*num_domain(2))

    ! Cell number of each dimension
    !
    cell(1) = boundary%num_cells_x
    cell(2) = boundary%num_cells_y
    cell(3) = boundary%num_cells_z

    ! Default value of the number of cell in each domain
    !
    ncell_local = 1

    ! Assign the cell index for each processor and the total number of cells
    !
    do i = 1, 3

      quotient = cell(i) / num_domain(i)
      remainder = mod(cell(i), num_domain(i))

      if (iproc(i) <= (remainder -1)) then
        quotient       = quotient + 1
        cell_start(i)  = quotient * iproc(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      else
        cell_start(i)  = (quotient+1)*remainder + quotient*(iproc(i)-remainder) 
        cell_start(i)  = cell_start(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      end if

    end do

    ! Assign the lower processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = num_domain(i) - 1

      iproc_lower(i) = iproc(1) + &
                       iproc(2)*num_domain(1) + &
                       iproc(3)*num_domain(1)*num_domain(2)

      iproc(i) = iproc(i) + 1
      if (iproc(i) == num_domain(i)) &
        iproc(i) = 0

    end do

    ! Assign the upper processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) + 1
      if (iproc(i) == num_domain(i)) &
        iproc(i) = 0

      iproc_upper(i) = iproc(1) + &
                       iproc(2)*num_domain(1) + &
                       iproc(3)*num_domain(1)*num_domain(2)

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = num_domain(i) - 1
    end do

    ! Assign the neighboring process index
    !
    do i = -1, 1

      if (i == -1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = num_domain(1) - 1

      else if (i == 1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == num_domain(1)) &
          iproc(1) = 0

      end if

      do j = -1, 1

        if (j == -1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = num_domain(2) - 1

        else if (j == 1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == num_domain(2)) &
            iproc(2) = 0

        end if

        do k = -1, 1

          if (k == -1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = num_domain(3) - 1

          else if (k == 1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == num_domain(3)) &
              iproc(3) = 0
          end if

          neighbor(i,j,k) = iproc(1) + &
                            iproc(2)*num_domain(1) + &
                            iproc(3)*num_domain(1)*num_domain(2)

          if (k == -1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == num_domain(3)) &
              iproc(3) = 0

          else if (k == 1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = num_domain(3) - 1

          end if
        end do

        if (j == -1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == num_domain(2)) &
            iproc(2) = 0

        else if (j == 1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = num_domain(2) - 1
        end if

      end do

      if (i == -1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == num_domain(1)) &
          iproc(1) = 0

      else if (i == 1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = num_domain(1) - 1

      end if
    end do

    ! maximum number of boundary cells
    !
    ncell_boundary = 2*(cell_length(1) * cell_length(2) +                   &
                        cell_length(2) * cell_length(3) +                   &
                        cell_length(1) * cell_length(3))+                   &
                     4*(cell_length(1) + cell_length(2) + cell_length(3)) + &
                     8
    return

  end subroutine setup_processor_rank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_capacity
  !> @brief        setup cell capacity for memory allocation
  !! @authors      JJ
  !! @param[in]    boundary : boundary information
  !! @param[in]    domain   : domain information
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_capacity(boundary, domain, molecule)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_molecule),        intent(in)    :: molecule

    ! local variables
    real(dp)                 :: v, v_rate
    integer                  :: natom, natom2, nbond, nangle
    integer                  :: ndihedral, nimproper, ncmap

    natom = molecule%num_atoms
    nbond = molecule%num_bonds
    nangle = molecule%num_angles
    ndihedral = molecule%num_dihedrals
    nimproper = molecule%num_impropers
    ncmap = molecule%num_cmaps

    v = boundary%box_size_x * &
        boundary%box_size_y * &
        boundary%box_size_z

    v = v / real(boundary%num_domain(1)* &
                 boundary%num_domain(2)* &
                 boundary%num_domain(3),dp)
    v = v / real(domain%num_cell_local,dp)

    v_rate = VboxRate*v / VolumeBox8

    MaxContact   = int(v_rate * real(NContBox8,wp) * ShrinkRate)
    ContactMove  = MaxContact / 2


    ! sp_domain_str
    !

    MaxAtom      = int(v_rate * real(NAtomBox8,dp))
#ifdef USE_GPU
    MaxAtom      = min(MaxAtom,NAtomMax_in_CUDA)
#endif
    MaxWater     = int(v_rate * real(NAtomBox8,dp) / 3.0_dp)
    MaxMove      = int(v_rate * real(NAtomBox8,dp) / 5.0_dp)
    MaxWaterMove = int(v_rate * real(NAtomBox8,dp) / 7.0_dp)

    ! If the system is vacuum, MaxAtom, MaxWater, etc. are set to 
    ! the number of atoms in the system.
    MaxAtom      = min(natom, MaxAtom)
    MaxWater     = min(natom, MaxWater)
    MaxMove      = min(natom, MaxMove)
    MaxWaterMove = min(natom, MaxWaterMove)


    ! sp_enefunc_str
    !

    MaxBond      = int(v_rate * real(NBondBox8,dp) * ShrinkRate)
    MaxAngle     = int(v_rate * real(NAnglBox8,dp) * ShrinkRate)
    MaxDihe      = int(v_rate * real(NDiheBox8,dp) * ShrinkRate)
    MaxImpr      = int(v_rate * real(NImprBox8,dp) * ShrinkRate)
    MaxCmap      = int(v_rate * real(NCmapBox8,dp) * ShrinkRate)

    BondMove     = MaxBond  / 2
    AngleMove    = MaxAngle / 2
    DiheMove     = MaxDihe  / 2
    ImprMove     = MaxImpr  / 2
    CmapMove     = MaxCmap  / 2

    ! If vacuum, MaxBond, MaxAngle, etc. are set to the number of
    ! bonds, angles, etc. of the molecule.
    MaxBond      = min(nbond, MaxBond)
    MaxAngle     = min(nangle, MaxAngle)
    MaxDihe      = min(10*ndihedral, MaxDihe)
    MaxImpr      = min(nimproper, MaxImpr)
    MaxCmap      = min(ncmap, MaxCmap)
    BondMove     = min(nbond, BondMove)
    AngleMove    = min(nangle, AngleMove)
    DiheMove     = min(10*ndihedral, DiheMove)
    ImprMove     = min(nimproper, ImprMove)
    CmapMove     = min(ncmap, CmapMove)


    ! sp_pairlist_str
    !

    MaxNb15      = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15      = min(natom, MaxNb15)
    MaxNb15      = MaxNb15 ** 2
    MaxNb15      = int(real(MaxNb15,dp) * Scale_pairlist)

    MaxNb15Water = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15Water = min(natom, MaxNb15Water)
    MaxNb15Water = MaxNb15Water ** 2
    MaxNb15water = MaxNb15Water / 10

    ! If vacuum, MaxNb15 and MaxNb15Water are set to natom**2.
    ! MaxContact and ContactMove are also set to natom**2.
    ! If the number of contact in the system is given, change natom**2 to it.
    MaxContact   = min(MaxNb15, MaxContact)
    ContactMove  = min(MaxNb15, ContactMove)


    ! sp_constraints_str
    !

    HGroupMax    = int(v_rate * real(NHGrpBox8,dp) * ShrinkRate)
    HGrpMaxMove  = int(v_rate * real(NHMovBox8,dp) * ShrinkRate)

    ! If vacuum, HGroupMax and HGrpMaxMove are set to natom.
    HGroupMax    = min(natom, HGroupMax)
    HGrpMaxMove  = min(natom, HGrpMaxMove)


#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) 'Cell volume                      : ', v
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxContact   : ', MaxContact
      write(MsgOut,*) 'sp_domain_str     ::ContactMove  : ', ContactMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxBond      : ', MaxBond
      write(MsgOut,*) 'sp_enefunc_str    ::MaxAngle     : ', MaxAngle
      write(MsgOut,*) 'sp_enefunc_str    ::MaxDihe      : ', MaxDihe
      write(MsgOut,*) 'sp_enefunc_str    ::MaxImpr      : ', MaxImpr
      write(MsgOut,*) 'sp_enefunc_str    ::MaxCmap      : ', MaxCmap
      write(MsgOut,*) 'sp_enefunc_str    ::BondMove     : ', BondMove
      write(MsgOut,*) 'sp_enefunc_str    ::AngleMove    : ', AngleMove
      write(MsgOut,*) 'sp_enefunc_str    ::DiheMove     : ', DiheMove
      write(MsgOut,*) 'sp_enefunc_str    ::ImprMove     : ', ImprMove
      write(MsgOut,*) 'sp_enefunc_str    ::CmapMove     : ', CmapMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_cell_capacity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_capacity_pio
  !> @brief        setup cell capacity for memory allocation for parallel io
  !! @authors      HO
  !! @param[in]    boundary : boundary information
  !! @param[in]    domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_capacity_pio(boundary, domain)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain

    ! local variables
    real(dp)                 :: v, v_rate
    
    v = boundary%box_size_x * &
        boundary%box_size_y * &
        boundary%box_size_z

    v = v / real(boundary%num_domain(1)* &
                 boundary%num_domain(2)* &
                 boundary%num_domain(3),dp)
    v = v / real(domain%num_cell_local,dp)

    v_rate = VboxRate*v / VolumeBox8

    MaxContact   = int(v_rate * real(NContBox8,wp) * ShrinkRate)
    ContactMove  = MaxContact / 2

    ! sp_domain_str
    !

    MaxAtom      = int(v_rate * real(NAtomBox8,dp))
#ifdef USE_GPU
    MaxAtom      = min(MaxAtom,NAtomMax_in_CUDA)
#endif
    MaxWater     = int(v_rate * real(NAtomBox8,dp) / 3.0_dp)
    MaxMove      = int(v_rate * real(NAtomBox8,dp) / 5.0_dp)
    MaxWaterMove = int(v_rate * real(NAtomBox8,dp) / 7.0_dp)


    ! sp_enefunc_str
    !

    MaxBond      = int(v_rate * real(NBondBox8,dp) * ShrinkRate)
    MaxAngle     = int(v_rate * real(NAnglBox8,dp) * ShrinkRate)
    MaxDihe      = int(v_rate * real(NDiheBox8,dp) * ShrinkRate)
    MaxImpr      = int(v_rate * real(NImprBox8,dp) * ShrinkRate)
    MaxCmap      = int(v_rate * real(NCmapBox8,dp) * ShrinkRate)

    BondMove     = MaxBond  / 2
    AngleMove    = MaxAngle / 2
    DiheMove     = MaxDihe  / 2
    ImprMove     = MaxImpr  / 2
    CmapMove     = MaxCmap  / 2


    ! sp_pairlist_str
    !

    MaxNb15      = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15      = MaxNb15 ** 2

    MaxNb15Water = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15Water = MaxNb15Water ** 2
    MaxNb15water = MaxNb15Water / 10


    ! sp_constraints_str
    !

    HGroupMax    = int(v_rate * real(NHGrpBox8,dp) * ShrinkRate)
    HGrpMaxMove  = int(v_rate * real(NHMovBox8,dp) * ShrinkRate)


#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) 'Cell volume                      : ', v
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxBond      : ', MaxBond
      write(MsgOut,*) 'sp_enefunc_str    ::MaxAngle     : ', MaxAngle
      write(MsgOut,*) 'sp_enefunc_str    ::MaxDihe      : ', MaxDihe
      write(MsgOut,*) 'sp_enefunc_str    ::MaxImpr      : ', MaxImpr
      write(MsgOut,*) 'sp_enefunc_str    ::MaxCmap      : ', MaxCmap
      write(MsgOut,*) 'sp_enefunc_str    ::BondMove     : ', BondMove
      write(MsgOut,*) 'sp_enefunc_str    ::AngleMove    : ', AngleMove
      write(MsgOut,*) 'sp_enefunc_str    ::DiheMove     : ', DiheMove
      write(MsgOut,*) 'sp_enefunc_str    ::ImprMove     : ', ImprMove
      write(MsgOut,*) 'sp_enefunc_str    ::CmapMove     : ', CmapMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_cell_capacity_pio


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_boundary
  !> @brief        define boundary cells in each domain
  !! @authors      JJ
  !! @param[in]    cell   : cell count for each axis
  !! @param[in]    ndom   : domain count for each axis
  !! @param[inout] domain : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_boundary(cell, ndom, domain)

    ! formal arguments
    integer,                 intent(in)    :: cell(3)
    integer,                 intent(in)    :: ndom(3)
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ic, jc, kc, icel, icel_local

    integer,         pointer :: ncell, ncell_boundary
    integer,         pointer :: cell_start(:), cell_end(:)
    integer,         pointer :: cell_g2b(:), cell_b2g(:)
    integer,         pointer :: cell_gxyz2l(:,:,:)
    integer,         pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,         pointer :: cell_l2gx_orig(:)
    integer,         pointer :: cell_l2gy_orig(:)
    integer,         pointer :: cell_l2gz_orig(:)


    ncell          => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end  
    cell_g2b       => domain%cell_g2b
    cell_b2g       => domain%cell_b2g
    cell_gxyz2l    => domain%cell_gxyz2l
    cell_l2gx      => domain%cell_l2gx
    cell_l2gy      => domain%cell_l2gy
    cell_l2gz      => domain%cell_l2gz
    cell_l2gx_orig => domain%cell_l2gx_orig
    cell_l2gy_orig => domain%cell_l2gy_orig
    cell_l2gz_orig => domain%cell_l2gz_orig

    icel_local = 0

    if (ndom(1) > 1) then
      ! face boundary (x direction, upper)
      !
      do j = cell_start(2), cell_end(2)
        do k = cell_start(3), cell_end(3)

          icel_local = icel_local + 1
          if (cell_end(1) == cell(1)) then
            i = 1
            ic = cell(1) + 1
          else
            i = cell_end(1) + 1
            ic = i
          end if

          icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = j
          cell_l2gz(icel_local+ncell) = k
          cell_l2gx_orig(icel_local+ncell) = i
          cell_l2gy_orig(icel_local+ncell) = j
          cell_l2gz_orig(icel_local+ncell) = k
          cell_gxyz2l(ic,j,k) = icel_local+ncell

        end do
      end do

      ! face boundary (x direction, lower)
      !
      do j = cell_start(2), cell_end(2)
        do k = cell_start(3), cell_end(3)

          icel_local = icel_local + 1
          if (cell_start(1) == 1) then
            i = cell(1)
            ic = 0
          else
            i = cell_start(1) - 1
            ic = i
          end if

          icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
          cell_g2b(icel) = icel_local
          cell_b2g(icel_local) = icel
          cell_l2gx(icel_local+ncell) = ic
          cell_l2gy(icel_local+ncell) = j
          cell_l2gz(icel_local+ncell) = k
          cell_l2gx_orig(icel_local+ncell) = i
          cell_l2gy_orig(icel_local+ncell) = j
          cell_l2gz_orig(icel_local+ncell) = k
          cell_gxyz2l(ic,j,k) = icel_local+ncell

        end do
      end do

      if (ndom(2) > 1) then
        ! face boundary (y direction, upper)
        !
        do ic = cell_start(1)-1, cell_end(1)+1

          if (ic == 0) then
            i = cell(1)
          else if (ic == (cell(1)+1)) then
            i = 1
          else
            i = ic
          end if

          do k = cell_start(3), cell_end(3)

            icel_local = icel_local + 1
            if (cell_end(2) == cell(2)) then
              j = 1
              jc = cell(2) + 1
            else
              j = cell_end(2) + 1
              jc = j
            end if

            icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
            cell_g2b(icel) = icel_local
            cell_b2g(icel_local) = icel
            cell_l2gx(icel_local+ncell) = ic
            cell_l2gy(icel_local+ncell) = jc
            cell_l2gz(icel_local+ncell) = k
            cell_l2gx_orig(icel_local+ncell) = i
            cell_l2gy_orig(icel_local+ncell) = j
            cell_l2gz_orig(icel_local+ncell) = k
            cell_gxyz2l(ic,jc,k) = icel_local+ncell

          end do
        end do

        ! face boundary (y direction, lower)
        !
        do ic = cell_start(1)-1, cell_end(1)+1

          if (ic == 0) then
            i = cell(1)
          else if (ic == (cell(1)+1)) then
            i = 1
          else
            i = ic
          end if

          do k = cell_start(3), cell_end(3)

            icel_local = icel_local + 1
            if (cell_start(2) == 1) then
              j = cell(2)
              jc = 0
            else
              j = cell_start(2) - 1
              jc = j
            end if

            icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
            cell_g2b(icel) = icel_local
            cell_b2g(icel_local) = icel
            cell_l2gx(icel_local+ncell) = ic
            cell_l2gy(icel_local+ncell) = jc
            cell_l2gz(icel_local+ncell) = k
            cell_l2gx_orig(icel_local+ncell) = i
            cell_l2gy_orig(icel_local+ncell) = j
            cell_l2gz_orig(icel_local+ncell) = k
            cell_gxyz2l(ic,jc,k) = icel_local+ncell

          end do
        end do

        if (ndom(3) > 1) then

          ! face boundary (z direction, upper)
          !
          do ic = cell_start(1)-1, cell_end(1)+1

            if (ic == 0) then
              i = cell(1)
            else if (ic == (cell(1)+1)) then
              i = 1
            else
              i = ic
            end if

            do jc = cell_start(2)-1, cell_end(2)+1

              if (jc == 0) then
                j = cell(2)
              else if (jc == (cell(2)+1)) then
                j = 1
              else
                j = jc
              end if

              icel_local = icel_local + 1
              if (cell_end(3) == cell(3)) then
                k = 1
                kc = cell(3) + 1
              else
                k = cell_end(3) + 1
                kc = k
              end if

              icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
              cell_g2b(icel) = icel_local
              cell_b2g(icel_local) = icel
              cell_l2gx(icel_local+ncell) = ic
              cell_l2gy(icel_local+ncell) = jc
              cell_l2gz(icel_local+ncell) = kc
              cell_l2gx_orig(icel_local+ncell) = i
              cell_l2gy_orig(icel_local+ncell) = j
              cell_l2gz_orig(icel_local+ncell) = k
              cell_gxyz2l(ic,jc,kc) = icel_local+ncell

            end do
          end do

          ! face boundary (z direction, lower)
          !
          do ic = cell_start(1)-1, cell_end(1)+1

            if (ic == 0) then
              i = cell(1)
            else if (ic == (cell(1)+1)) then
              i = 1
            else
              i = ic
            end if

            do jc = cell_start(2)-1, cell_end(2)+1

              if (jc == 0) then
                j = cell(2)
              else if (jc == (cell(2)+1)) then
                j = 1
              else
                j = jc
              end if

              icel_local = icel_local + 1
              if (cell_start(3) == 1) then
                k = cell(3)
                kc = 0
              else
                k = cell_start(3) - 1
                kc = k
              end if

              icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
              cell_g2b(icel) = icel_local
              cell_b2g(icel_local) = icel
              cell_l2gx(icel_local+ncell) = ic
              cell_l2gy(icel_local+ncell) = jc
              cell_l2gz(icel_local+ncell) = kc
              cell_l2gx_orig(icel_local+ncell) = i
              cell_l2gy_orig(icel_local+ncell) = j
              cell_l2gz_orig(icel_local+ncell) = k
              cell_gxyz2l(ic,jc,kc) = icel_local+ncell

            end do
          end do
        end if
      end if
    end if

    ! total number of boundary cells
    !
    ncell_boundary = icel_local

    return

  end subroutine setup_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_nobc_boundary
  !> @brief        update system size in NOBC condition
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[inout] boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_nobc_boundary(domain, boundary)

    ! formal arguments
    type(s_domain),           intent(in)    :: domain
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    integer         :: i, ix
    real(dp)        :: coord_min(1:3), coord_max(1:3), box_size(1:3)

    coord_min(1:3)          =  1000000000000.0_dp
    coord_max(1:3)          = -1000000000000.0_dp
    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        coord_min(1:3) = min(coord_min(1:3), domain%coord(1:3,ix,i))
        coord_max(1:3) = max(coord_max(1:3), domain%coord(1:3,ix,i))
      end do
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, coord_min, 3, mpi_real8, mpi_min, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, coord_max, 3, mpi_real8, mpi_max, &
                       mpi_comm_country, ierror)
#endif
    box_size(1:3) = max(-coord_min(1:3), coord_max(1:3)) + 0.1_dp
    boundary%box_size_x = box_size(1)*2.0_dp
    boundary%box_size_y = box_size(2)*2.0_dp
    boundary%box_size_z = box_size(3)*2.0_dp
    if (boundary%box_size_x > boundary%box_size_x_max) &
      call error_msg('calculated box_size_x is greater than box_size_x_max')
    if (boundary%box_size_y > boundary%box_size_y_max) &
      call error_msg('calculated box_size_y is greater than box_size_y_max')
    if (boundary%box_size_z > boundary%box_size_z_max) &
      call error_msg('calculated box_size_z is greater than box_size_z_max')
    if (boundary%box_size_x < boundary%box_size_x_min) &
      write(Msgout, '(A)') 'WARNING : calculated box size_x is less than box_size_x_min'
    if (boundary%box_size_y < boundary%box_size_y_min) &
      write(Msgout, '(A)') 'WARNING : calculated box size_y is less than box_size_y_min'
    if (boundary%box_size_z < boundary%box_size_z_min) &
      write(Msgout, '(A)') 'WARNING : calculated box size_z is less than box_size_z_min'

    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,wp)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,wp)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,wp)

    return

  end subroutine update_nobc_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_and_water
  !> @brief        setup solute and water list
  !! @authors      JJ
  !! @param[in]    water_model : water model
  !! @param[in]    molecule    : molecule information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_and_water(molecule, enefunc, water_model, tip4)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    character(*),            intent(inout) :: water_model
    logical,                 intent(inout) :: tip4

    ! local variables
    integer                  :: i, k, natom, nwater, nsolute, water_atom
    integer                  :: io, ih, id(4)
    real(wp)                 :: mass(4), mass_max, mass_min


    natom   = molecule%num_atoms
    nwater  = 0
    nsolute = 0

    do i = 1, natom
      if ((molecule%residue_name(i)(1:3) == 'TIP' .or.  &
           molecule%residue_name(i)(1:3) == 'WAT' .or.  &
           molecule%residue_name(i)(1:3) == 'SOL' .or.  &
           molecule%residue_name(i)(1:3) == 'SPC').and. &
           .not.(molecule%light_atom_mass(i))     .and. &
           .not.(molecule%light_atom_name(i))) then
        nwater = nwater + 1
      else if (molecule%residue_name(i)(1:3) /= 'TIP' .and. &
               molecule%residue_name(i)(1:3) /= 'WAT' .and. &
               molecule%residue_name(i)(1:3) /= 'SOL' .and. &
               molecule%residue_name(i)(1:3) /= 'SPC') then
        nsolute = nsolute + 1
      end if
    end do

    if ((natom-nsolute)/4 == nwater .and. nwater > 0) then
      tip4 = .true.
    else if ((natom-nsolute)/3 == nwater) then
      tip4 = .false.
    else
      call error_msg('Setup_Solute_And_Water> # of water is incorrect.')
    end if

    enefunc%table%num_water  = nwater
    enefunc%table%num_solute = nsolute
    call alloc_enefunc(enefunc, EneFuncTableWat, nwater,  1)
    call alloc_enefunc(enefunc, EneFuncTableSol, nsolute, natom)

    if (tip4) then
      water_atom = 4
    else
      water_atom = 3
    end if

    i = 1
    nwater  = 0
    nsolute = 0

    do while(.true.)

      if (i > natom) exit

      if (molecule%residue_name(i)(1:3) == 'TIP' .or. &
          molecule%residue_name(i)(1:3) == 'WAT' .or. &
          molecule%residue_name(i)(1:3) == 'SOL' .or. &
          molecule%residue_name(i)(1:3) == 'SPC') then

        do k = 1, water_atom
          mass(k) = molecule%mass(i-1+k)
        end do

        mass_max = -1000.0_wp
        mass_min = 1000.0_wp
        do k = 1, water_atom
          mass_max = max(mass_max,mass(k))
          mass_min = min(mass_min, mass(k))
        end do

        id(1:water_atom) = 0

        if (water_atom == 4) then
          do k = 1, water_atom
            if (mass(k) == mass_max) id(1) = i-1+k
            if (mass(k) == mass_min) id(water_atom) = i-1+k
            if (mass(k) > mass_min .and. mass(k) < mass_max) then
              if (id(2) == 0) id(2) = i-1+k
              if (id(2) /= 0) id(3) = i-1+k
            end if
          end do
        else if (water_atom == 3) then
          do k = 1, water_atom
            if (mass(k) == mass_max) id(1) = i-1+k
            if (mass(k) == mass_min) then
              if (id(2) == 0) id(2) = i-1+k
              if (id(2) /= 0) id(3) = i-1+k
            end if
          end do
        end if

        nwater = nwater + 1
        enefunc%table%water_list(1:water_atom,nwater) = id(1:water_atom)
        i = i + water_atom

      else

        if (molecule%mass(i) .lt. EPS) &
          call error_msg('Setup_Solute_And_Water> Mass = 0 is not allowed')

        nsolute = nsolute + 1
        enefunc%table%solute_list(nsolute) = i
        enefunc%table%solute_list_inv(i) = nsolute
        i = i + 1

      end if

    end do

    if (nwater /= enefunc%table%num_water) &
      call error_msg('Setup_Solute_And_Water> number of water is incorrect')

    if (nsolute /= enefunc%table%num_solute) &
      call error_msg('Setup_Solute_And_Water> number of solute is incorrect')


    ! set water oxygen and hydrogen
    !
    if (size(enefunc%table%water_list(1,:)) /= 0) then

      io = enefunc%table%water_list(1,1)
      ih = enefunc%table%water_list(2,1)

      enefunc%table%atom_cls_no_O = molecule%atom_cls_no(io)
      enefunc%table%atom_cls_no_H = molecule%atom_cls_no(ih)
      enefunc%table%charge_O      = molecule%charge(io)
      enefunc%table%charge_H      = molecule%charge(ih)
      enefunc%table%mass_O        = molecule%mass(io)
      enefunc%table%mass_H        = molecule%mass(ih)

      ! dummy
      !
      if (tip4) then
        io = enefunc%table%water_list(4,1)
        enefunc%table%atom_cls_no_D = molecule%atom_cls_no(io)
        enefunc%table%charge_D      = molecule%charge(io)
        enefunc%table%mass_D        = molecule%mass(io)
      end if

    end if

    return

  end subroutine setup_solute_and_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hbond_group
  !> @brief        setup bonds including hydrogen atom
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hbond_group(molecule, enefunc, constraints)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, l, i1, i2, nhmax, nbond_p
    character(6)             :: ci1, ci2
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2

    ! count the number of bonds including hydrogen
    !
    call alloc_constraints(constraints, ConstraintsHBond, molecule%num_atoms)

    do i = 1, molecule%num_bonds

      i1 = molecule%bond_list(1,i)
      i2 = molecule%bond_list(2,i)
      mi1 = molecule%light_atom_mass(i1)
      mi2 = molecule%light_atom_mass(i2)
      cl1 = molecule%light_atom_name(i1)
      cl2 = molecule%light_atom_name(i2)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      endif

      if (enefunc%table%solute_list_inv(i1) /= 0 .and. &
          enefunc%table%solute_list_inv(i2) /= 0) then

!        ci1 = molecule%atom_name(i1)
!        ci2 = molecule%atom_name(i2)

        if (cl1 .or. cl2) then

          if (cl1) then
            constraints%duplicate(i2) = constraints%duplicate(i2) + 1
            constraints%H_index(constraints%duplicate(i2),i2) = i1
          else
            constraints%duplicate(i1) = constraints%duplicate(i1) + 1
            constraints%H_index(constraints%duplicate(i1),i1) = i2
          end if

        end if

      end if

    end do


    ! count XHn group for each number n
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then
        constraints%nh(1) = constraints%nh(1) + 1

      else if (constraints%duplicate(i1) == 2) then
        constraints%nh(2) = constraints%nh(2) + 1

      else if (constraints%duplicate(i1) == 3) then
        constraints%nh(3) = constraints%nh(3) + 1

      else if (constraints%duplicate(i1) == 4) then
        constraints%nh(4) = constraints%nh(4) + 1

      else if (constraints%duplicate(i1) == 5) then
        constraints%nh(5) = constraints%nh(5) + 1

      else if (constraints%duplicate(i1) == 6) then
        constraints%nh(6) = constraints%nh(6) + 1

      else if (constraints%duplicate(i1) == 7) then
        constraints%nh(7) = constraints%nh(7) + 1

      else if (constraints%duplicate(i1) == 8) then
        constraints%nh(8) = constraints%nh(8) + 1

      else if (constraints%duplicate(i1) >= 8) then
        call error_msg( &
             'Setup_HBond_Group> Bond(>8) for one atom is not considered')

      end if

    end do

    constraints%connect  = 0
    if (constraints%nh(1) /= 0) constraints%connect = 1
    if (constraints%nh(2) /= 0) constraints%connect = 2
    if (constraints%nh(3) /= 0) constraints%connect = 3
    if (constraints%nh(4) /= 0) constraints%connect = 4
    if (constraints%nh(5) /= 0) constraints%connect = 5
    if (constraints%nh(6) /= 0) constraints%connect = 6
    if (constraints%nh(7) /= 0) constraints%connect = 7
    if (constraints%nh(8) /= 0) constraints%connect = 8

    nhmax = max(constraints%nh(1),constraints%nh(2),constraints%nh(3), &
                constraints%nh(4),constraints%nh(5),constraints%nh(6), &
                constraints%nh(7),constraints%nh(8))

    call alloc_constraints(constraints, ConstraintsBondGroup, nhmax)


    ! Make a list of XHn
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then

        constraints%nh(1) = constraints%nh(1) + 1
        constraints%H_Group(1,constraints%nh(1),1) = i1
        constraints%H_Group(2,constraints%nh(1),1) = constraints%H_index(1,i1)

      else if (constraints%duplicate(i1) == 2) then

        constraints%nh(2) = constraints%nh(2) + 1
        constraints%H_Group(1,constraints%nh(2),2) = i1
        constraints%H_Group(2,constraints%nh(2),2) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(2),2) = constraints%H_index(2,i1)

      else if (constraints%duplicate(i1) == 3) then

        constraints%nh(3) = constraints%nh(3) + 1
        constraints%H_Group(1,constraints%nh(3),3) = i1
        constraints%H_Group(2,constraints%nh(3),3) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(3),3) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(3),3) = constraints%H_index(3,i1)

      else if (constraints%duplicate(i1) == 4) then

        constraints%nh(4) = constraints%nh(4) + 1
        constraints%H_Group(1,constraints%nh(4),4) = i1
        constraints%H_Group(2,constraints%nh(4),4) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(4),4) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(4),4) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(4),4) = constraints%H_index(4,i1)

      else if (constraints%duplicate(i1) == 5) then

        constraints%nh(5) = constraints%nh(5) + 1
        constraints%H_Group(1,constraints%nh(5),5) = i1
        constraints%H_Group(2,constraints%nh(5),5) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(5),5) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(5),5) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(5),5) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(5),5) = constraints%H_index(5,i1)

      else if (constraints%duplicate(i1) == 6) then

        constraints%nh(6) = constraints%nh(6) + 1
        constraints%H_Group(1,constraints%nh(6),6) = i1
        constraints%H_Group(2,constraints%nh(6),6) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(6),6) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(6),6) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(6),6) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(6),6) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(6),6) = constraints%H_index(6,i1)

      else if (constraints%duplicate(i1) == 7) then

        constraints%nh(7) = constraints%nh(7) + 1
        constraints%H_Group(1,constraints%nh(7),7) = i1
        constraints%H_Group(2,constraints%nh(7),7) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(7),7) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(7),7) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(7),7) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(7),7) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(7),7) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(7),7) = constraints%H_index(7,i1)

      else if (constraints%duplicate(i1) == 8) then

        constraints%nh(8) = constraints%nh(8) + 1
        constraints%H_Group(1,constraints%nh(8),8) = i1
        constraints%H_Group(2,constraints%nh(8),8) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(8),8) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(8),8) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(8),8) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(8),8) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(8),8) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(8),8) = constraints%H_index(7,i1)
        constraints%H_Group(9,constraints%nh(8),8) = constraints%H_index(8,i1)

      end if

    end do


    return

  end subroutine setup_hbond_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_HBond
  !> @brief        setup atom maps with H-bond connection groups
  !! @authors      JJ
  !! @param[in]    molecule    : molecule information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_HBond(molecule, boundary, enefunc, constraints, &
                                 domain)

    ! formal arguments
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_enefunc),     target, intent(in)    :: enefunc
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(dp)                     :: x_shift, y_shift, z_shift
    real(dp)                     :: move(3), origin(3)
    integer                      :: i, j, icx, icy, icz, icel
    integer                      :: isolute, iwater, ih1, ih2, id
    integer                      :: icel_local
    integer                      :: ncel_local, ncel
    integer                      :: patm, psol, pwat, pnoh, phgl
    character(4)                 :: ci1
    logical                      :: mi1, cl1, tip4

    real(dp),            pointer :: bsize_x, bsize_y, bsize_z
    real(dp),            pointer :: csize_x, csize_y, csize_z
    integer,             pointer :: ncel_x, ncel_y, ncel_z
    integer,             pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: natom(:), nsolute(:), nwater(:)
    integer,             pointer :: solute_list(:,:), water_list(:,:,:)
    integer,             pointer :: No_HGr(:), HGr_local(:,:)
    integer,             pointer :: HGr_bond_list(:,:,:,:)

    ncel          = domain%num_cell_local + domain%num_cell_boundary
    ncel_local    = domain%num_cell_local

    call alloc_constraints(constraints, ConstraintsDomainBond, ncel, &
                           constraints%connect)

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water
    solute_list   => domain%solute_list
    water_list    => domain%water_list

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list

    tip4          = constraints%tip4
    origin(1)     = boundary%origin_x
    origin(2)     = boundary%origin_y
    origin(3)     = boundary%origin_z


    ! solute atoms (not bonded to hydrogen) in each domain
    !
    do isolute = 1, enefunc%table%num_solute

      i   = enefunc%table%solute_list(isolute)
      ci1 = molecule%atom_name(i)

      mi1 = molecule%light_atom_mass(i)
      cl1 = molecule%light_atom_name(i)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
      endif

!      if (constraints%duplicate(i) == 0 .and. ci1(1:1) /= 'H') then
      if (constraints%duplicate(i) == 0 .and. .not. cl1) then

        !coordinate shifted against the origin
        !
        x_shift = molecule%atom_coord(1,i) - boundary%origin_x
        y_shift = molecule%atom_coord(2,i) - boundary%origin_y
        z_shift = molecule%atom_coord(3,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        ! atoms inside the domain
        !
        if (cell_g2l(icel) /= 0) then

          ! local cell index
          !
          icel_local = cell_g2l(icel)

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          pnoh =  No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          natom  (icel_local) = patm
          nsolute(icel_local) = psol
          No_HGr (icel_local) = pnoh

        ! atoms in the boundary
        !
        else if (cell_g2b(icel) /= 0) then
  
          ! local cell index
          !
          icel_local = cell_g2b(icel) + ncel_local

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          pnoh =  No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          natom  (icel_local) = patm
          nsolute(icel_local) = psol
          No_HGr (icel_local) = pnoh

        end if 

      end if

    end do

    ! Hydrogen bonding group
    !
    do isolute = 1, constraints%connect

      do j = 1, constraints%nh(isolute)

        i = constraints%H_Group(1,j,isolute)
        x_shift = molecule%atom_coord(1,i) - boundary%origin_x
        y_shift = molecule%atom_coord(2,i) - boundary%origin_y
        z_shift = molecule%atom_coord(3,i) - boundary%origin_z

        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        ! atoms inside the domain
        !
        if (cell_g2l(icel) /= 0) then
  
          icel_local = cell_g2l(icel)

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          phgl =  HGr_local(isolute,icel_local)
  
          ! hydrogen atom
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm

          call molecule_to_domain(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          HGr_local(isolute,icel_local) = phgl

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1,j,isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm

            call molecule_to_domain(molecule, move, origin, i, &
                                    domain, icel_local, patm)

          end do

          natom  (icel_local) = patm
          nsolute(icel_local) = psol

        else if (cell_g2b(icel) /= 0) then

          icel_local = cell_g2b(icel) + ncel_local

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          phgl =  HGr_local(isolute,icel_local)

          ! hydrogen atoms
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm

          call molecule_to_domain(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          HGr_local(isolute,icel_local) = phgl

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1,j,isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm

            call molecule_to_domain(molecule, move, origin, i, &
                                    domain, icel_local, patm)

          end do

          natom  (icel_local) = patm
          nsolute(icel_local) = psol

        end if

      end do

    end do

    ! water atoms in each domain
    !
    do iwater = 1, enefunc%table%num_water

      i   = enefunc%table%water_list(1,iwater)
      ih1 = enefunc%table%water_list(2,iwater)
      ih2 = enefunc%table%water_list(3,iwater)
      if (tip4) id = enefunc%table%water_list(4,iwater)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm
          call molecule_to_domain(molecule, move, origin, id, &
                                  domain, icel_local, patm)
        end if 
          
        natom(icel_local)  = patm
        nwater(icel_local) = pwat

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm
          call molecule_to_domain(molecule, move, origin, id, &
                                  domain, icel_local, patm)
        end if

        natom(icel_local)  = patm
        nwater(icel_local) = pwat

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_HBond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_table
  !> @brief        setup atom maps with solute and water table
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_table(molecule, boundary, enefunc, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3), origin(3)
    integer                   :: i, icx, icy, icz, icel
    integer                   :: isolute, iwater, ih1, ih2, id
    integer                   :: icel_local
    integer                   :: ncel_local, ncel
    integer                   :: patm, psol, pwat
    logical                   :: tip4

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: solute_list(:,:), water_list(:,:,:)


    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    ncel_x       => boundary%num_cells_x
    ncel_y       => boundary%num_cells_y
    ncel_z       => boundary%num_cells_z
    csize_x      => boundary%cell_size_x
    csize_y      => boundary%cell_size_y
    csize_z      => boundary%cell_size_z

    cell_g2l     => domain%cell_g2l
    cell_g2b     => domain%cell_g2b
    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water
    solute_list  => domain%solute_list
    water_list   => domain%water_list

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z

    tip4         = enefunc%table%tip4

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local


    ! solute atoms in each domain
    !
    do isolute = 1, enefunc%table%num_solute

      i = enefunc%table%solute_list(isolute)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)
      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom  (icel_local)
        psol =  nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1
        solute_list(psol,icel_local) = patm
        solute_list(patm,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        natom  (icel_local) = patm
        nsolute(icel_local) = psol

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom  (icel_local)
        psol =  nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1
        solute_list(psol,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        natom  (icel_local) = patm
        nsolute(icel_local) = psol

      end if

    end do

    ! water atoms in each domain
    !
    do iwater = 1, enefunc%table%num_water

      i   = enefunc%table%water_list(1,iwater)
      ih1 = enefunc%table%water_list(2,iwater)
      ih2 = enefunc%table%water_list(3,iwater)
      if (tip4) id = enefunc%table%water_list(4,iwater)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm
          call molecule_to_domain(molecule, move, origin, id, &
                                  domain, icel_local, patm)
        end if

        natom (icel_local) = patm
        nwater(icel_local) = pwat

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm
          call molecule_to_domain(molecule, move, origin, id, &
                                  domain, icel_local, patm)
        end if

        natom (icel_local) = patm
        nwater(icel_local) = pwat

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom
  !> @brief        setup atom maps with whole atoms
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom(molecule, boundary, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3), origin(3)
    integer                   :: i, icx, icy, icz, icel
    integer                   :: icel_local
    integer                   :: ncel_local, ncel, natom_all

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), nsolute(:)


    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    csize_x      => boundary%cell_size_x
    csize_y      => boundary%cell_size_y
    csize_z      => boundary%cell_size_z
    ncel_x       => boundary%num_cells_x
    ncel_y       => boundary%num_cells_y
    ncel_z       => boundary%num_cells_z

    cell_g2l     => domain%cell_g2l
    cell_g2b     => domain%cell_g2b
    natom        => domain%num_atom
    nsolute      => domain%num_solute

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local
    natom_all    = domain%num_atom_all


    do i = 1, natom_all

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        natom(icel_local) = natom(icel_local) + 1

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, natom(icel_local))

      ! atoms in a boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        natom(icel_local) = natom(icel_local) + 1

        call molecule_to_domain(molecule, move, origin, i, &
                                domain, icel_local, natom(icel_local))

      end if

    end do

    nsolute(1:ncel) = natom(1:ncel)
    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_neighboring_cells
  !> @brief        check the neighboring cells of each cell
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_neighbor_cells(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, j, ii, ic, jc, inbc, k
    integer                   :: ncel_local, nboundary

    integer,          pointer :: cell_g2l(:), cell_l2g(:)
    integer,          pointer :: cell_b2g(:), cell_g2b(:)


    cell_g2l   => domain%cell_g2l
    cell_g2b   => domain%cell_g2b
    cell_b2g   => domain%cell_b2g
    cell_l2g   => domain%cell_l2g

    ncel_local = domain%num_cell_local
    nboundary  = domain%num_cell_boundary

    call alloc_domain(domain, DomainNeighbourCell, ncel_local+nboundary, 1, 1)

    do i = 1, ncel_local

      k  = 0
      ic = cell_l2g(i)

      do inbc = 1, 27

        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2l(jc)

        if (j > i) then

          k = k + 1
          domain%neighbor_cells(k,i) = j

        else if (cell_g2b(jc) /= 0) then

          j = cell_g2b(jc) + ncel_local
          k = k + 1
          domain%neighbor_cells(k,i) = j

        end if

      end do

      domain%near_neighbor_cell_count(i) = k

      do inbc = 28, 125

        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2l(jc)

        if (j > i) then

          k = k + 1
          domain%neighbor_cells(k,i) = j

        else if (cell_g2b(jc) /= 0) then

          j = cell_g2b(jc) + ncel_local
          k = k + 1
          domain%neighbor_cells(k,i) = j

        end if

      end do

      domain%neighbor_cell_count(i) = k

    end do

    do ii = 1, nboundary

      k  = 0
      i  = ii + ncel_local
      ic = cell_b2g(ii)

      do inbc = 1, 27

        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2b(jc) + ncel_local

        if (j > i) then

          k = k + 1
          domain%neighbor_cells(k,i) = j

        end if
      end do

      domain%near_neighbor_cell_count(i) = k

      do inbc = 28, 125

        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2b(jc) + ncel_local 

        if (j > i) then

          k = k + 1
          domain%neighbor_cells(k,i) = j

        end if
      end do

      domain%neighbor_cell_count(i) = k

    end do

    return

  end subroutine assign_neighbor_cells
    
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_atoms
  !> @brief        compate the total number of atoms in each cell
  !! @authors      JJ
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_atoms(natom, num_atom,                  &
                               cell_l2gx,  cell_l2gy, cell_l2gz, &
                               cell_start, cell_end,  neighbor,  &
                               ncel_local, nboundary)

    ! formal arguments
    real(wp),                intent(inout) :: natom(:)
    integer,                 intent(in)    :: num_atom(:) 
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: neighbor(-1:1,-1:1,-1:1)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary

    ! local variables
    integer                  :: i, ii, i1, i2, i3 


    do i = 1, ncel_local
      natom(i) = real(num_atom(i)*nproc_city + my_city_rank, wp)
    end do

    do ii = 1, nboundary

      i  = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)

      if (i1 == cell_start(1)-1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,0),wp)
          end if
        end if
      else if (i1 == cell_end(1)+1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,0),wp)
          end if
        end if
      else if (i1 >= cell_start(1) .and. i1 <= cell_end(1)) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,0),wp)
          end if
        end if
      end if
    end do

    return

  end subroutine assign_cell_atoms


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_cpu  
  !> @brief        compute the total cput time of each cell (randomized)
  !! @authors      JJ
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_cpu(time_proc, random,               &
                             cell_l2gx, cell_l2gy, cell_l2gz, &
                             cell_start, cell_end, neighbor,  &
                             ncel_local, nboundary, cpu_time)

    ! formal arguments
    real(wp),                intent(in)    :: time_proc(:)
    real(dp),                intent(in)    :: random(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: neighbor(-1:1,-1:1,-1:1)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary
    real(wp),                intent(inout) :: cpu_time(:)

    ! local variables
    integer                  :: i, ii, i1, i2, i3, j1, j2, j3


    do i = 1, ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)
      j1 = i1 + 1
      j2 = i2 + 1
      j3 = i3 + 1
      cpu_time(i) = time_proc(my_city_rank+1)*random(i)
    end do

    do ii = 1, nboundary

      i = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)
      j1 = i1 + 1
      j2 = i2 + 1
      j3 = i3 + 1

      if (i1 == cell_start(1)-1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,0,0)+1)*random(i)
          end if
        end if
      else if (i1 == cell_end(1)+1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,0,0)+1)*random(i)
          end if
        end if
      else if (i1 >= cell_start(1) .and. i1 <= cell_end(1)) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,0,0)+1)*random(i)
          end if
        end if
      end if
    end do

    return

  end subroutine assign_cell_cpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_interactions
  !> @brief        assign the cell index for given cell-cell interaction
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_interaction(natom, bsize_x, bsize_y, bsize_z,   &
                                     num_domain, cell,                   &
                                     cell_l2gx,  cell_l2gy,  cell_l2gz,  &
                                     cell_l2gx1, cell_l2gy1, cell_l2gz1, &
                                     cell_gxyz2l,cell_start, cell_end,   &
                                     ncel_local, nboundary,              &
                                     cell_move, cell_pair, virial_check)

    ! formal arguments
    real(wp),                intent(in)    :: natom(:)
    real(dp),                intent(in)    :: bsize_x
    real(dp),                intent(in)    :: bsize_y
    real(dp),                intent(in)    :: bsize_z
    integer,                 intent(in)    :: num_domain(:)
    integer,                 intent(in)    :: cell(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_l2gx1(:)
    integer,                 intent(in)    :: cell_l2gy1(:)
    integer,                 intent(in)    :: cell_l2gz1(:)
    integer,                 intent(in)    :: cell_gxyz2l(:,:,:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary
    real(wp),                intent(inout) :: cell_move(:,:,:)
    integer,                 intent(inout) :: cell_pair(:,:)
    integer,                 intent(inout) :: virial_check(:,:)

    ! local variables
    real(dp)                 :: ic1, ic2
    integer                  :: i, i1, i2, i3, j, ii, jj, ij
    integer                  :: icx1, icy1, icz1, icx2, icy2, icz2
    integer                  :: icx, icy, icz, movex, movey, movez


    ij = 0

    ! assign the interaction cell for each interaction
    !
    do i = 1, ncel_local-1

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do j = i+1, ncel_local

        cell_pair(i,j) = 0
        cell_pair(j,i) = 0

        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then
          icx = (icx1 + icx2)/2
          icy = (icy1 + icy2)/2
          icz = (icz1 + icz2)/2
          cell_pair(i,j) = i
          cell_pair(j,i) = i

          icx =min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
          icy =min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
          icz =min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

          if (icx == abs(icx1-icx2-cell(1))) cell_move(1,j,i) = -1.0_wp
          if (icx == abs(icx1-icx2+cell(1))) cell_move(1,j,i) =  1.0_wp
          if (icx == abs(icx1-icx2))         cell_move(1,j,i) =  0.0_wp
          if (icy == abs(icy1-icy2-cell(2))) cell_move(2,j,i) = -1.0_wp
          if (icy == abs(icy1-icy2+cell(2))) cell_move(2,j,i) =  1.0_wp
          if (icy == abs(icy1-icy2))         cell_move(2,j,i) =  0.0_wp
          if (icz == abs(icz1-icz2-cell(3))) cell_move(3,j,i) = -1.0_wp
          if (icz == abs(icz1-icz2+cell(3))) cell_move(3,j,i) =  1.0_wp
          if (icz == abs(icz1-icz2))         cell_move(3,j,i) =  0.0_wp
          if (abs(cell_move(1,j,i)) <= EPS .and. &
              abs(cell_move(2,j,i)) <= EPS .and. &
              abs(cell_move(3,j,i)) <= EPS) then
            virial_check(j,i) = 0
            virial_check(i,j) = 0
          end if

          cell_move(1:3,i,j) = - cell_move(1:3,j,i)

          ij = ij + 1
        end if

      end do
    end do

    do i = 1, ncel_local
      cell_pair(i,i) = i
      ij = ij + 1
    end do

    do i = 1, ncel_local

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do jj = 1, nboundary

        j = jj + ncel_local
        cell_pair(j,i) = 0
        cell_pair(i,j) = 0

        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        movex = 0
        movey = 0
        movez = 0

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (icx == abs(icx1-icx2-cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = -1.0_wp
          movex = -cell(1)
        else if (icx == abs(icx1-icx2+cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = 1.0_wp
          movex = cell(1)
        else
          cell_move(1,j,i) = 0.0_wp
          movex = 0
        end if

        if (icy == abs(icy1-icy2-cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = -1.0_wp
          movey = -cell(2)
        else if (icy == abs(icy1-icy2+cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = 1.0_wp
          movey = cell(2)
        else
          cell_move(2,j,i) = 0.0_wp
          movey = 0
        end if

        if (icz == abs(icz1-icz2-cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = -1.0_wp
          movez = -cell(3)
        else if (icz == abs(icz1-icz2+cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = 1.0_wp
          movez = cell(3)
        else
          cell_move(3,j,i) = 0.0_wp
          movez = 0
        end if

        cell_move(1:3,i,j) = - cell_move(1:3,j,i)

        icx = icx1 - icx2 + movex
        icy = icy1 - icy2 + movey
        icz = icz1 - icz2 + movez

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

          ij = ij + 1
          icx = icx1 + icx2 + movex
          icy = icy1 + icy2 + movey
          icz = icz1 + icz2 + movez

          if (icx == (2*cell_start(1)-1) .or. icx == (2*cell_end(1)+1)) then

            ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
            ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))

            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else
            i1 = icx / 2
          end if

          if ((icy == (2*cell_start(2)-1) .or. icy == (2*cell_end(2)+1))) then
            if (num_domain(2) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1)) 
              if (ic1 < ic2) then
                i2 = icy1
              else
                i2 = icy2
              end if
            else
              i2 = (icy1 + icy2) / 2
            end if
          else
            i2 = (icy1 + icy2) / 2
          end if

          if ((icz == (2*cell_start(3)-1) .or. icz == (2*cell_end(3)+1))) then
            if (num_domain(3) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i3 = icz1
              else
                i3 = icz2
              end if
            else
              i3 = (icz1 + icz2) / 2
            end if
          else
            i3 = (icz1 + icz2) / 2
          end if

          cell_pair(i,j) = cell_gxyz2l(i1+1,i2+1,i3+1)
          cell_pair(j,i) = cell_gxyz2l(i1+1,i2+1,i3+1)

        end if

      end do
    end do

    do i = 1, ncel_local

      icx1 = cell_l2gx1(i)
      icy1 = cell_l2gy1(i)
      icz1 = cell_l2gz1(i)

      do jj = 1, nboundary

        j = jj + ncel_local

        icx2 = cell_l2gx1(j)
        icy2 = cell_l2gy1(j)
        icz2 = cell_l2gz1(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

        if (icx == abs(icx1-icx2-cell(1))) then
          cell_move(1,j,i) = -1.0_wp
        else if (icx == abs(icx1-icx2+cell(1))) then
          cell_move(1,j,i) = 1.0_wp
        else
          cell_move(1,j,i) = 0.0_wp
        end if

        if (icy == abs(icy1-icy2-cell(2))) then
          cell_move(2,j,i) = -1.0_wp
        else if (icy == abs(icy1-icy2+cell(2))) then
          cell_move(2,j,i) = 1.0_wp
        else
          cell_move(2,j,i) = 0.0_wp
        end if

        if (icz == abs(icz1-icz2-cell(3))) then
          cell_move(3,j,i) = -1.0_wp
        else if (icz == abs(icz1-icz2+cell(3))) then
          cell_move(3,j,i) = 1.0_wp
        else
          cell_move(3,j,i) = 0.0_wp
        end if
        cell_move(1:3,i,j) = -cell_move(1:3,j,i)

        if (abs(cell_move(1,j,i)) <= EPS .and. &
            abs(cell_move(2,j,i)) <= EPS .and. &
            abs(cell_move(3,j,i)) <= EPS) then
          virial_check(j,i) = 0
          virial_check(i,j) = 0
        end if

        end if

      end do
    end do

    do ii = 1, nboundary
      ij = ij + 1
      i = ii + ncel_local
      cell_pair(i,i) = i
    end do

    do ii = 1, nboundary-1

      i = ii + ncel_local
      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do jj = ii+1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)
        cell_pair(j,i) = 0
        cell_pair(i,j) = 0

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        movex = 0
        movey = 0
        movez = 0
        icx = min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (icx == abs(icx1-icx2-cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = -1.0_wp
          movex = -cell(1)
        else if (icx == abs(icx1-icx2+cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = 1.0_wp
          movex = cell(1)
        else
          cell_move(1,j,i) = 0.0_wp
          movex = 0
        end if

        if (icy == abs(icy1-icy2-cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = -1.0_wp
          movey = -cell(2)
        else if (icy == abs(icy1-icy2+cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = 1.0_wp
          movey = cell(2)
        else
          cell_move(2,j,i) = 0.0_wp
          movey = 0
        end if

        if (icz == abs(icz1-icz2-cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = -1.0_wp
          movez = -cell(3)
        else if (icz == abs(icz1-icz2+cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = 1.0_wp
          movez = cell(3)
        else
          cell_move(3,j,i) = 0.0_wp
          movez = 0
        end if

        cell_move(1:3,i,j) = - cell_move(1:3,j,i)

        icx = icx1 - icx2 + movex
        icy = icy1 - icy2 + movey
        icz = icz1 - icz2 + movez

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

          ij = ij + 1
          icx = icx1 + icx2
          icy = icy1 + icy2 + movey
          icz = icz1 + icz2 + movez

          if (icx == (2*cell_start(1)-1) .or. icx == (2*cell_end(1)+1)) then
            ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
            ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else
            i1 = icx / 2
          end if

          if ((icy == (2*cell_start(2)-1) .or. icy == (2*cell_end(2)+1))) then
            if (num_domain(2) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i2 = icy1
              else
                i2 = icy2
              end if
            else
              i2 = (icy1 + icy2) / 2
            end if
          else
            i2 = (icy1 + icy2) / 2
          end if

          if ((icz == (2*cell_start(3)-1) .or. icz == (2*cell_end(3)+1))) then
            if (num_domain(3) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i3 = icz1
              else
                i3 = icz2
              end if
            else
              i3 = (icz1 + icz2) / 2
            end if
          else
            i3 = (icz1 + icz2) / 2
          end if

          cell_pair(i,j) = cell_gxyz2l(i1+1,i2+1,i3+1)
          cell_pair(j,i) = cell_gxyz2l(i1+1,i2+1,i3+1)

        end if
      end do
    end do

    do ii = 1, nboundary-1

      i = ii + ncel_local
      icx1 = cell_l2gx1(i)
      icy1 = cell_l2gy1(i)
      icz1 = cell_l2gz1(i)

      do jj = ii+1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx1(j)
        icy2 = cell_l2gy1(j)
        icz2 = cell_l2gz1(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

        if (icx == abs(icx1-icx2-cell(1))) then
          cell_move(1,j,i) = -1.0_wp
        else if (icx == abs(icx1-icx2+cell(1))) then
          cell_move(1,j,i) = 1.0_wp
        else
          cell_move(1,j,i) = 0.0_wp
        end if

        if (icy == abs(icy1-icy2-cell(2))) then
          cell_move(2,j,i) = -1.0_wp
        else if (icy == abs(icy1-icy2+cell(2))) then
          cell_move(2,j,i) = 1.0_wp
        else
          cell_move(2,j,i) = 0.0_wp
        end if

        if (icz == abs(icz1-icz2-cell(3))) then
          cell_move(3,j,i) = -1.0_wp
        else if (icz == abs(icz1-icz2+cell(3))) then
          cell_move(3,j,i) = 1.0_wp
        else
          cell_move(3,j,i) = 0.0_wp
        end if

        cell_move(1,i,j) = -cell_move(1,j,i)
        cell_move(2,i,j) = -cell_move(2,j,i)
        cell_move(3,i,j) = -cell_move(3,j,i)

        if (abs(cell_move(1,j,i)) <= EPS .and. &
            abs(cell_move(2,j,i)) <= EPS .and. &
            abs(cell_move(3,j,i)) <= EPS) then
          virial_check(j,i) = 0
          virial_check(i,j) = 0
        end if

        end if

      end do
    end do

    univ_maxcell1 = ij + ncel_local

    return
  
  end subroutine assign_cell_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain
  !> @brief        copy molecule information to domain
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain(molecule, move, origin, iatom, &
                                domain, icel, icel_atom)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    real(dp),                 intent(in)    :: move(3)
    real(dp),                 intent(in)    :: origin(3)
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    ! local variables
    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)

    real(dp),         pointer :: coord_local (:,:,:)
    real(dp),         pointer :: ref_local   (:,:,:)
    real(dp),         pointer :: vel_local   (:,:,:)
    real(wp),         pointer :: charge_local(:,:)
    real(dp),         pointer :: mass_local  (:,:)
    real(wp),         pointer :: trans       (:,:,:)
    integer,          pointer :: class_local (:,:)
    integer,          pointer :: id_l2g      (:,:)
    integer,          pointer :: id_g2l      (:,:)
    

    coord        => molecule%atom_coord
    vel          => molecule%atom_velocity
    charge       => molecule%charge
    mass         => molecule%mass
    atom_class   => molecule%atom_cls_no

    id_g2l       => domain%id_g2l
    id_l2g       => domain%id_l2g
    coord_local  => domain%coord
    ref_local    => domain%coord_ref
    vel_local    => domain%velocity
    charge_local => domain%charge
    mass_local   => domain%mass
    class_local  => domain%atom_cls_no
    trans        => domain%trans_vec


    id_l2g(icel_atom,icel) = iatom
    id_g2l(1,iatom) = icel
    id_g2l(2,iatom) = icel_atom


    coord_local(1:3,icel_atom,icel) = coord (1:3,iatom) - origin(1:3)
    ref_local  (1:3,icel_atom,icel) = coord (1:3,iatom) - origin(1:3)
    vel_local  (1:3,icel_atom,icel) = vel   (1:3,iatom)
    charge_local   (icel_atom,icel) = charge    (iatom)
    mass_local     (icel_atom,icel) = mass      (iatom)
    class_local    (icel_atom,icel) = atom_class(iatom)
    trans      (1:3,icel_atom,icel) = move(1:3) 

    return

  end subroutine molecule_to_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord
  !> @brief        check_atom_coordinate
  !! @authors      CK
  !! @param[in]    ene_info      : ENERGY section control parameters information
  !! @param[in]    boundary      : boundary condition information
  !! @param[in]    contact_check : flag for contact_check
  !! @param[inout] domain        : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord(ene_info, boundary, contact_check, domain)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    type(s_boundary), target, intent(in)    :: boundary
    logical,                  intent(in)    :: contact_check
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    integer                   :: i, ix, iy, ij, j
    integer                   :: id, my_id,  omp_get_thread_num
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: cell_pairlist1(:,:)
    integer,          pointer :: id_l2g(:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(dp),         pointer :: mass(:,:), coord(:,:,:) 
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)

    mass           => domain%mass 
    coord          => domain%coord
    trans1         => domain%trans_vec
    trans2         => domain%translated
    cell_move      => domain%cell_move
    system_size    => domain%system_size
    natom          => domain%num_atom
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    cell_pairlist1 => domain%cell_pairlist1
    id_l2g         => domain%id_l2g

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
        end do
      end do
    else
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) 
        end do
      end do
    endif
    !$omp barrier

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        dir(1:3) = trans2(1:3,ix,i)
        do iy = ix + 1, natom(i)
          dr(1:3) = dir(1:3) - trans2(1:3,iy,i)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          if (rr < EPS) then
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
            call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
          end if

          if (rr < ene_info%minimum_contact .and. &
              abs(mass(ix,i)) > EPS .and. abs(mass(iy,i)) > EPS) then 
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
            if (rr < ene_info%err_minimum_contact .and.  &
               .not. contact_check) &
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
          endif
        end do
      end do
    end do

    if (boundary%type == BoundaryTypePBC) then
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)

        tr(1:3) = cell_move(1:3,j,i)*system_size(1:3)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i) 
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - coord(1:3,iy,j)+ tr(1:3)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < ene_info%minimum_contact .and. &
                abs(mass(ix,i)) > EPS .and. abs(mass(iy,j)) > EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and.  &
                 .not. contact_check) &
                call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
            endif
     
          end do
        end do
      end do
    else
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i)
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - trans2(1:3,iy,j)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < ene_info%minimum_contact .and. &
                abs(mass(ix,i)) > EPS .and. abs(mass(iy,j)) > EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and.  &
                 .not. contact_check) &
                call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
            endif
     
          end do
        end do
      end do
    endif

    !$omp end parallel

    return

  end subroutine check_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord_pio
  !> @brief        check_atom_coordinate
  !! @authors      CK
  !! @param[in]    enefunc  : energy potential function information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord_pio(enefunc, boundary, domain)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    integer                   :: i, ix, iy, ij, j
    integer                   :: id, my_id,  omp_get_thread_num
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: cell_pairlist1(:,:)
    integer,          pointer :: id_l2g(:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(dp),         pointer :: coord(:,:,:) 
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)

    coord          => domain%coord
    trans1         => domain%trans_vec
    trans2         => domain%translated
    cell_move      => domain%cell_move
    system_size    => domain%system_size
    natom          => domain%num_atom
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    cell_pairlist1 => domain%cell_pairlist1
    id_l2g         => domain%id_l2g

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
        end do
      end do
    else
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) 
        end do
      end do
    endif
    !$omp barrier

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        dir(1:3) = trans2(1:3,ix,i)
        do iy = ix + 1, natom(i)
          dr(1:3) = dir(1:3) - trans2(1:3,iy,i)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          if (rr < EPS) then
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
            call error_msg('Check_Atom_Coord_Pio> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
          end if

          if (rr < enefunc%minimum_contact) then 
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
          endif
        end do
      end do
    end do

    if (boundary%type == BoundaryTypePBC) then
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)

        tr(1:3) = cell_move(1:3,j,i)*system_size(1:3)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i) 
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - coord(1:3,iy,j)+ tr(1:3)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord_Pio> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < enefunc%minimum_contact) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
            endif
     
          end do
        end do
      end do
    else
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i)
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - trans2(1:3,iy,j)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord_Pio> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < enefunc%minimum_contact) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
            endif
     
          end do
        end do
      end do
    endif

    !$omp end parallel

    return

  end subroutine check_atom_coord_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_solute_and_water
  !> @brief        compate the total number of atoms in each cell
  !! @authors      JJ
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_solute_and_water(natom, nsolute, nwater,          &
                                          cell_l2gx, cell_l2gy, cell_l2gz, &
                                          cell_start, cell_end,            &
                                          neighbor, ncel_local, nboundary)

    ! formal arguments
    real(wp),                intent(inout) :: natom(:)
    integer,                 intent(in)    :: nsolute(:)
    integer,                 intent(in)    :: nwater(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: neighbor(-1:1,-1:1,-1:1)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary

    ! local variables
    integer                  :: i, ii, i1, i2, i3


    do i = 1, ncel_local
      natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
      natom(i) = natom(i)*real(nproc_city,wp) + real(my_city_rank,wp)
    end do

    do ii = 1, nboundary

      i = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)

      if (i1 == cell_start(1)-1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(-1,0,0),wp)
          end if
        end if
      else if (i1 == cell_end(1)+1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(1,0,0),wp)
          end if
        end if
      else if (i1 >= cell_start(1) .and. i1 <= cell_end(1)) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(nsolute(i),wp) + 1.0_wp*real(nwater(i),wp)
            natom(i) = natom(i)*real(nproc_city,wp)+real(neighbor(0,0,0),wp)
          end if
        end if
      end if

    end do

    return

  end subroutine assign_cell_solute_and_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_fep
  !> @brief        setup domain information for FEP calculation
  !! @authors      NK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_fep(ene_info, con_info, &
                          boundary, molecule, enefunc, constraints, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all


    ! initialize structure informations
    !
    call init_domain(domain)
    call init_enefunc(enefunc)
    call init_constraints(constraints)

    domain%num_atom_all       = molecule%num_atoms

    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type

    domain%num_atom_single_all  = molecule%num_atoms_fep(2)
    domain%fep_use     = .true.

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)


    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity_fep(boundary, domain, molecule)


    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)


    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do


    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)


    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar_FEP, ncel_all, 1, 1)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all,    1, 1)

    ! decide hydrogen atom from mass
    !
    call check_light_atom_name(con_info%hydrogen_mass_upper_bound, molecule)

    if (boundary%type == BoundaryTypePBC) then

      if (constraints%rigid_bond) then

        call setup_solute_and_water(molecule, enefunc,       &
                                    constraints%water_model, &
                                    constraints%tip4)

        if (constraints%tip4) then
          call alloc_domain(domain, DomainDynvar_Atom_FEP, ncel_all, 4, 1)
        else
          call alloc_domain(domain, DomainDynvar_Atom_FEP, ncel_all, 3, 1)
        end if

        call setup_hbond_group_fep     (molecule, enefunc, constraints)

        call setup_atom_by_HBond_fep   (molecule, boundary, enefunc, &
                                    constraints, domain)

      else 

        call setup_solute_and_water(molecule, enefunc,         &
                                    enefunc%table%water_model, &
                                    enefunc%table%tip4)

        if (enefunc%table%tip4) then
          call alloc_domain(domain, DomainDynvar_Atom_FEP, ncel_all, 4, 1)
        else
          call alloc_domain(domain, DomainDynvar_Atom_FEP, ncel_all, 3, 1)
        end if
        call setup_atom_by_table_fep (molecule, boundary, enefunc, domain)

      end if

    else if (boundary%type == BoundaryTypeNOBC) then

      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 0, 1)
      call setup_atom_fep (molecule, boundary, domain)

    end if

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    call check_atom_coord_fep(ene_info, boundary, ene_info%contact_check, domain)

    ! sort atoms assigned to single topology
    call sort_single_fep(domain, enefunc)

    ! setup water molecule information
    !
    if (boundary%type == BoundaryTypePBC .and. &
        enefunc%table%num_water > 0) then

      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
      domain%water%charge(1)      = enefunc%table%charge_O
      domain%water%charge(2)      = enefunc%table%charge_H
      domain%water%charge(3)      = enefunc%table%charge_H
      domain%water%charge(4)      = enefunc%table%charge_D
      domain%water%mass(1)        = enefunc%table%mass_O
      domain%water%mass(2)        = enefunc%table%mass_H
      domain%water%mass(3)        = enefunc%table%mass_H
      domain%water%mass(4)        = enefunc%table%mass_D

    end if

    return

  end subroutine setup_domain_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_ptl_fep
  !> @brief        update particles in each cell for FEP calculations
  !! @authors      NK
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_ptl_fep(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3)
    integer                   :: i, k, ix, icx, icy, icz, icel, ncel
    integer                   :: icel_local, icel_bd

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    real(dp),         pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(dp),         pointer :: mass(:,:)
    real(dp),         pointer :: buf_real(:,:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: atmcls(:,:), id_l2g(:,:)
    integer,          pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,          pointer :: buf_int(:,:,:)
    integer,          pointer :: fepgrp(:,:)


    bsize_x  => boundary%box_size_x
    bsize_y  => boundary%box_size_y
    bsize_z  => boundary%box_size_z
    ncel_x   => boundary%num_cells_x
    ncel_y   => boundary%num_cells_y
    ncel_z   => boundary%num_cells_z
    csize_x  => boundary%cell_size_x
    csize_y  => boundary%cell_size_y
    csize_z  => boundary%cell_size_z

    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real

    ! FEP
    fepgrp   => domain%fepgrp

    ! initializaiton
    !
    ncel = domain%num_cell_local + domain%num_cell_boundary
    ptl_exit(1:domain%num_cell_local) = 0
    ptl_add (1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, domain%num_cell_local

      k = 0
      do ix = 1, domain%num_atom(i)

        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1

        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        icel_local = domain%cell_g2l(icel)
        icel_bd    = domain%cell_g2b(icel)

        if (icel_local /= i) then

          ptl_exit(i) = ptl_exit(i) + 1
          ptlindex(ptl_exit(i),i) = ix

          if (icel_local /= 0) then

            ptl_add(icel_local) = ptl_add(icel_local) + 1
            buf_real(1,ptl_add(icel_local),icel_local) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_local),icel_local) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_local),icel_local) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_local),icel_local) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_local),icel_local) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_local),icel_local) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_local),icel_local) = charge(ix,i)
            buf_real(8,ptl_add(icel_local),icel_local) = mass(ix,i)
            buf_int (1,ptl_add(icel_local),icel_local) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_local),icel_local) = id_l2g(ix,i)
            buf_int (3,ptl_add(icel_local),icel_local) = fepgrp(ix,i)

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + domain%num_cell_local 
            ptl_add(icel_bd) = ptl_add(icel_bd) + 1
            buf_real(1,ptl_add(icel_bd),icel_bd) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_bd),icel_bd) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_bd),icel_bd) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_bd),icel_bd) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_bd),icel_bd) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_bd),icel_bd) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_bd),icel_bd) = charge(ix,i)
            buf_real(8,ptl_add(icel_bd),icel_bd) = mass(ix,i)
            buf_int (1,ptl_add(icel_bd),icel_bd) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_bd),icel_bd) = id_l2g(ix,i)
            buf_int (3,ptl_add(icel_bd),icel_bd) = fepgrp(ix,i)

          end if
        end if
      end do
    end do

    return

  end subroutine update_outgoing_ptl_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_ptl_fep
  !> @brief        update particles in each cell for FEP
  !! @authors      NK
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_ptl_fep(domain)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(dp),        pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: buf_real(:,:,:)
    integer,         pointer :: atmcls(:,:), id_l2g(:,:), id_g2l(:,:)
    integer,         pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,         pointer :: buf_int(:,:,:)
    integer,         pointer :: fepgrp(:,:)


    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    id_g2l   => domain%id_g2l
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real

    ! FEP
    fepgrp   => domain%fepgrp


    ! Incoming particles
    !
    do i = 1, domain%num_cell_local

      ! When the number of coming particles is larger than that of outgoing ones
      !
      if (ptl_add(i) >= ptl_exit(i)) then

        do k = 1, ptl_exit(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          fepgrp(ptlindex(k,i),i)     = buf_int(3,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        do k = ptl_exit(i)+1, ptl_add(i)
          ix = k + domain%num_atom(i) - ptl_exit(i)
          coord(1,ix,i)               = buf_real(1,k,i)
          coord(2,ix,i)               = buf_real(2,k,i)
          coord(3,ix,i)               = buf_real(3,k,i)
          velocity(1,ix,i)            = buf_real(4,k,i)
          velocity(2,ix,i)            = buf_real(5,k,i)
          velocity(3,ix,i)            = buf_real(6,k,i)
          charge(ix,i)                = buf_real(7,k,i)
          mass(ix,i)                  = buf_real(8,k,i)
          atmcls(ix,i)                = buf_int(1,k,i)
          id_l2g(ix,i)                = buf_int(2,k,i)
          fepgrp(ix,i)                = buf_int(3,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ix
        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)

      ! When the number of coming particles is less than that of outgoing ones
      !
      else

        do k = 1, ptl_add(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          fepgrp(ptlindex(k,i),i)     = buf_int(3,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        j  = 0
        ix = domain%num_atom(i)
        k  = ptl_add(i) + 1

        do while (j < (ptl_exit(i)-ptl_add(i)))

          insert = .true.
          do kx = k, ptl_exit(i)
            if (ix == ptlindex(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = ptlindex(k,i)
            coord(1,kx,i)          = coord(1,ix,i)
            coord(2,kx,i)          = coord(2,ix,i)
            coord(3,kx,i)          = coord(3,ix,i)
            velocity(1,kx,i)       = velocity(1,ix,i)
            velocity(2,kx,i)       = velocity(2,ix,i)
            velocity(3,kx,i)       = velocity(3,ix,i)
            charge(kx,i)           = charge(ix,i)
            mass(kx,i)             = mass(ix,i)
            atmcls(kx,i)           = atmcls(ix,i)
            id_l2g(kx,i)           = id_l2g(ix,i)
            fepgrp(kx,i)           = fepgrp(ix,i)
            id_g2l(1,id_l2g(kx,i)) = i
            id_g2l(2,id_l2g(kx,i)) = kx

            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)
           
      end if

      domain%num_solute(i) = domain%num_atom(i)

    end do

    return

  end subroutine update_incoming_ptl_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_capacity_fep
  !> @brief        setup cell capacity for memory allocation for FEP
  !! @authors      HO
  !! @param[in]    boundary : boundary information
  !! @param[in]    domain   : domain information
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_capacity_fep(boundary, domain, molecule)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_molecule),        intent(in)    :: molecule

    ! local variables
    real(dp)                 :: v, v_rate
    integer                  :: natom, natom2, nbond, nangle
    integer                  :: ndihedral, nimproper, ncmap

    ! FEP
    integer                  :: natom_partA
    integer                  :: natom_partB
  
    natom = molecule%num_atoms
    nbond = molecule%num_bonds
    nangle = molecule%num_angles
    ndihedral = molecule%num_dihedrals
    nimproper = molecule%num_impropers
    ncmap = molecule%num_cmaps

    v = boundary%box_size_x * &
        boundary%box_size_y * &
        boundary%box_size_z

    v = v / real(boundary%num_domain(1)* &
                 boundary%num_domain(2)* &
                 boundary%num_domain(3),dp)
    v = v / real(domain%num_cell_local,dp)

    v_rate = VboxRate*v / VolumeBox8

    ! FEP
    ! Since in FEP perturbed atoms can overlap water molecules,
    ! each degree of freedom rarely exceeds the "Max***".
    ! If not sufficient, change this ratio.
    v_rate = v_rate * 1.4_dp

    MaxContact   = int(v_rate * real(NContBox8,wp) * ShrinkRate)
    ContactMove  = MaxContact / 2


    ! sp_domain_str
    !

    MaxAtom      = int(v_rate * real(NAtomBox8,dp))
#ifdef USE_GPU
    MaxAtom      = min(MaxAtom,NAtomMax_in_CUDA)
#endif
    MaxWater     = int(v_rate * real(NAtomBox8,dp) / 3.0_dp)
    MaxMove      = int(v_rate * real(NAtomBox8,dp) / 5.0_dp)
    MaxWaterMove = int(v_rate * real(NAtomBox8,dp) / 7.0_dp)

    ! If the system is vacuum, MaxAtom, MaxWater, etc. are set to 
    ! the number of atoms in the system.
    MaxAtom      = min(natom, MaxAtom)
    MaxWater     = min(natom, MaxWater)
    MaxMove      = min(natom, MaxMove)
    MaxWaterMove = min(natom, MaxWaterMove)


    ! sp_enefunc_str
    !

    MaxBond      = int(v_rate * real(NBondBox8,dp) * ShrinkRate)
    MaxAngle     = int(v_rate * real(NAnglBox8,dp) * ShrinkRate)
    MaxDihe      = int(v_rate * real(NDiheBox8,dp) * ShrinkRate)
    MaxImpr      = int(v_rate * real(NImprBox8,dp) * ShrinkRate)
    MaxCmap      = int(v_rate * real(NCmapBox8,dp) * ShrinkRate)

    BondMove     = MaxBond  / 2
    AngleMove    = MaxAngle / 2
    DiheMove     = MaxDihe  / 2
    ImprMove     = MaxImpr  / 2
    CmapMove     = MaxCmap  / 2

    ! If vacuum, MaxBond, MaxAngle, etc. are set to the number of
    ! bonds, angles, etc. of the molecule.
    MaxBond      = min(nbond, MaxBond)
    MaxAngle     = min(nangle, MaxAngle)
    MaxDihe      = min(10*ndihedral, MaxDihe)
    MaxImpr      = min(nimproper, MaxImpr)
    MaxCmap      = min(ncmap, MaxCmap)
    BondMove     = min(nbond, BondMove)
    AngleMove    = min(nangle, AngleMove)
    DiheMove     = min(10*ndihedral, DiheMove)
    ImprMove     = min(nimproper, ImprMove)
    CmapMove     = min(ncmap, CmapMove)


    ! sp_pairlist_str
    !

    MaxNb15      = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15      = min(natom, MaxNb15)

    ! FEP
    ! Number of pairs between perturbed atoms and others
    ! is, at most, MaxNb15 * (Number of perturbed atoms).
    natom_partA = molecule%num_atoms_fep(1)+molecule%num_atoms_fep(2)
    natom_partB = molecule%num_atoms_fep(3)+molecule%num_atoms_fep(4)
    MaxNb15_fep = max(natom_partA, natom_partB)
    MaxNb15_fep = min(MaxNb15, MaxNb15_fep)
    MaxNb15_fep = MaxNb15 * MaxNb15_fep

    MaxNb15      = MaxNb15 ** 2

    MaxNb15Water = int(v_rate * real(NAtomBox8,dp) * ShrinkRate)
    MaxNb15Water = min(natom, MaxNb15Water)
    MaxNb15Water = MaxNb15Water ** 2
    MaxNb15water = MaxNb15Water / 10

    ! If vacuum, MaxNb15 and MaxNb15Water are set to natom**2.
    ! MaxContact and ContactMove are also set to natom**2.
    ! If the number of contact in the system is given, change natom**2 to it.
    MaxContact   = min(MaxNb15, MaxContact)
    ContactMove  = min(MaxNb15, ContactMove)


    ! sp_constraints_str
    !

    HGroupMax    = int(v_rate * real(NHGrpBox8,dp) * ShrinkRate)
    HGrpMaxMove  = int(v_rate * real(NHMovBox8,dp) * ShrinkRate)

    ! If vacuum, HGroupMax and HGrpMaxMove are set to natom.
    HGroupMax    = min(natom, HGroupMax)
    HGrpMaxMove  = min(natom, HGrpMaxMove)


#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) 'Cell volume                      : ', v
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxContact   : ', MaxContact
      write(MsgOut,*) 'sp_domain_str     ::ContactMove  : ', ContactMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxBond      : ', MaxBond
      write(MsgOut,*) 'sp_enefunc_str    ::MaxAngle     : ', MaxAngle
      write(MsgOut,*) 'sp_enefunc_str    ::MaxDihe      : ', MaxDihe
      write(MsgOut,*) 'sp_enefunc_str    ::MaxImpr      : ', MaxImpr
      write(MsgOut,*) 'sp_enefunc_str    ::MaxCmap      : ', MaxCmap
      write(MsgOut,*) 'sp_enefunc_str    ::BondMove     : ', BondMove
      write(MsgOut,*) 'sp_enefunc_str    ::AngleMove    : ', AngleMove
      write(MsgOut,*) 'sp_enefunc_str    ::DiheMove     : ', DiheMove
      write(MsgOut,*) 'sp_enefunc_str    ::ImprMove     : ', ImprMove
      write(MsgOut,*) 'sp_enefunc_str    ::CmapMove     : ', CmapMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_fep  : ', MaxNb15_fep
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_cell_capacity_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pmelist_fep
  !> @brief        update pmelist in each cell for FEP calculations
  !! @authors      NK
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pmelist_fep(domain, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc

    ! local variable
    integer                   :: i, ix

    integer,          pointer :: natom(:)
    integer,          pointer :: natom_preserve(:)
    integer,          pointer :: natom_singleA(:)
    integer,          pointer :: natom_singleB(:)
    integer,          pointer :: natom_vanish_gr(:)
    integer,          pointer :: natom_appear_gr(:)
    integer,          pointer :: pmelist_preserve(:,:)
    integer,          pointer :: pmelist_appear_gr(:,:)
    integer,          pointer :: pmelist_vanish_gr(:,:)
    integer,          pointer :: id_singleA(:,:,:)
    integer,          pointer :: id_singleB(:,:,:)
    integer                   :: fg


    natom             => domain%num_atom
    natom_preserve    => domain%num_atom_preserve
    natom_singleA     => domain%num_atom_singleA
    natom_singleB     => domain%num_atom_singleB
    natom_appear_gr   => domain%num_atom_appear_gr
    natom_vanish_gr   => domain%num_atom_vanish_gr
    pmelist_preserve  => domain%pmelist_preserve
    pmelist_appear_gr => domain%pmelist_appear_gr
    pmelist_vanish_gr => domain%pmelist_vanish_gr
    id_singleA        => domain%id_singleA
    id_singleB        => domain%id_singleB

    ! Initialize
    !
    natom_preserve(:)      = 0
    natom_singleA(:)       = 0
    natom_singleB(:)       = 0
    natom_appear_gr(:)     = 0
    natom_vanish_gr(:)     = 0
    pmelist_preserve(:,:)  = 0
    pmelist_appear_gr(:,:) = 0
    pmelist_vanish_gr(:,:) = 0
    id_singleA(:,:,:)      = 0
    id_singleB(:,:,:)      = 0

    ! Check outgoing particles
    !
    do i = 1, domain%num_cell_local + domain%num_cell_boundary

      do ix = 1, natom(i)

        fg = domain%fepgrp(ix,i)
        if(enefunc%fep_topology == 2) then
          if (fg == 5) then
            natom_preserve(i)  = natom_preserve(i)  + 1
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_preserve(natom_preserve(i),i)   = ix
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 1) then
            natom_preserve(i)  = natom_preserve(i)  + 1
            natom_singleA(i)   = natom_singleA(i)   + 1
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_preserve(natom_preserve(i),i)   = ix
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 2) then
            natom_preserve(i)  = natom_preserve(i)  + 1
            natom_singleB(i)   = natom_singleB(i)   + 1
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_preserve(natom_preserve(i),i)   = ix
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 3) then
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 4) then
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
          end if
        else
          if (fg == 5) then
            natom_preserve(i)  = natom_preserve(i)  + 1
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_preserve(natom_preserve(i),i)   = ix
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 1) then
            natom_singleA(i)   = natom_singleA(i)   + 1
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
            id_singleA(natom_singleA(i),i,1) = ix
            id_singleA(natom_singleA(i),i,2) = domain%id_l2g(ix, i)
            id_singleA(natom_singleA(i),i,3) = natom_singleA(i)
          else if (fg == 2) then
            natom_singleB(i)   = natom_singleB(i)   + 1
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
            id_singleB(natom_singleB(i),i,1) = ix
            id_singleB(natom_singleB(i),i,2) = domain%id_l2g(ix, i)
            id_singleB(natom_singleB(i),i,3) = natom_singleB(i)
          else if (fg == 3) then
            natom_vanish_gr(i) = natom_vanish_gr(i) + 1
            pmelist_vanish_gr(natom_vanish_gr(i),i) = ix
          else if (fg == 4) then
            natom_appear_gr(i) = natom_appear_gr(i) + 1
            pmelist_appear_gr(natom_appear_gr(i),i) = ix
          end if
        end if

      end do

    end do

    return

  end subroutine update_pmelist_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hbond_group_fep
  !> @brief        setup bonds including hydrogen atom for FEP
  !! @authors      HO
  !! @param[in]    molecule : molecule information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hbond_group_fep(molecule, enefunc, constraints)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, l, i1, i2, nhmax, nbond_p
    character(6)             :: ci1, ci2
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2

    ! count the number of bonds including hydrogen
    !
    call alloc_constraints(constraints, ConstraintsHBond, molecule%num_atoms)

    do i = 1, molecule%num_bonds

      i1 = molecule%bond_list(1,i)
      i2 = molecule%bond_list(2,i)
      mi1 = molecule%light_atom_mass(i1)
      mi2 = molecule%light_atom_mass(i2)
      cl1 = molecule%light_atom_name(i1)
      cl2 = molecule%light_atom_name(i2)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1 
        cl2 = mi2 
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
        cl2 = (cl2 .or. mi2) 
      endif

      if (enefunc%table%solute_list_inv(i1) /= 0 .and. &
          enefunc%table%solute_list_inv(i2) /= 0) then

!        ci1 = molecule%atom_name(i1)
!        ci2 = molecule%atom_name(i2)

        if (cl1 .or. cl2) then

          ! Hydrogen rewiring in FEP:
          ! When singleA-dualA or singleB-dualB bonds have hydrogen atoms,
          ! SHAKE problems occur. If coordinate and velocity are synchonized
          ! after SHAKE, coordinate and velocity of atoms involved in the bonds
          ! are moved and then the constraints are sometimes broken. To avoid
          ! this problem, hydrogen atoms beloging to dualB is rewired to connect
          ! with atoms of singleA. The constraints of the hydrogen atoms are
          ! considered in performing SHAKE for molecule A. After SHAKE, singleB
          ! atoms are synchonized to singleA.
          !
          ! singleA--dualA
          !                 ==>   singleA--dualA
          ! singleB--dualB               \_dualB
          ! 
          ! To rewire the bonds, for singleA-dualB bond including hydrogen,
          ! the atom index of singleB is replaced with the corresponding atom
          ! index of singleA.
          if ((int(molecule%fepgrp(i1)) == 2) .and. &
              (int(molecule%fepgrp(i2)) == 4)) then
            do k = 1, molecule%num_atoms_fep(1)
              if (molecule%id_singleB(k) == i1) then
                i1 = molecule%id_singleA(k)
                exit
              end if
            end do
          else if ((int(molecule%fepgrp(i1)) == 4) .and. &
                   (int(molecule%fepgrp(i2)) == 2)) then
            do k = 1, molecule%num_atoms_fep(2)
              if (molecule%id_singleB(k) == i2) then
                i2 = molecule%id_singleA(k)
                exit
              end if
            end do
          end if

          if (cl1) then
            constraints%duplicate(i2) = constraints%duplicate(i2) + 1
            constraints%H_index(constraints%duplicate(i2),i2) = i1
          else
            constraints%duplicate(i1) = constraints%duplicate(i1) + 1
            constraints%H_index(constraints%duplicate(i1),i1) = i2
          end if

        end if

      end if

    end do


    ! count XHn group for each number n
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then
        constraints%nh(1) = constraints%nh(1) + 1

      else if (constraints%duplicate(i1) == 2) then
        constraints%nh(2) = constraints%nh(2) + 1

      else if (constraints%duplicate(i1) == 3) then
        constraints%nh(3) = constraints%nh(3) + 1

      else if (constraints%duplicate(i1) == 4) then
        constraints%nh(4) = constraints%nh(4) + 1

      else if (constraints%duplicate(i1) == 5) then
        constraints%nh(5) = constraints%nh(5) + 1

      else if (constraints%duplicate(i1) == 6) then
        constraints%nh(6) = constraints%nh(6) + 1

      else if (constraints%duplicate(i1) == 7) then
        constraints%nh(7) = constraints%nh(7) + 1

      else if (constraints%duplicate(i1) == 8) then
        constraints%nh(8) = constraints%nh(8) + 1

      else if (constraints%duplicate(i1) >= 8) then
        call error_msg( &
             'Setup_HBond_Group> Bond(>8) for one atom is not considered')

      end if

    end do

    constraints%connect  = 0
    if (constraints%nh(1) /= 0) constraints%connect = 1
    if (constraints%nh(2) /= 0) constraints%connect = 2
    if (constraints%nh(3) /= 0) constraints%connect = 3
    if (constraints%nh(4) /= 0) constraints%connect = 4
    if (constraints%nh(5) /= 0) constraints%connect = 5
    if (constraints%nh(6) /= 0) constraints%connect = 6
    if (constraints%nh(7) /= 0) constraints%connect = 7
    if (constraints%nh(8) /= 0) constraints%connect = 8

    nhmax = max(constraints%nh(1),constraints%nh(2),constraints%nh(3), &
                constraints%nh(4),constraints%nh(5),constraints%nh(6), &
                constraints%nh(7),constraints%nh(8))

    call alloc_constraints(constraints, ConstraintsBondGroup, nhmax)


    ! Make a list of XHn
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then

        constraints%nh(1) = constraints%nh(1) + 1
        constraints%H_Group(1,constraints%nh(1),1) = i1
        constraints%H_Group(2,constraints%nh(1),1) = constraints%H_index(1,i1)

      else if (constraints%duplicate(i1) == 2) then

        constraints%nh(2) = constraints%nh(2) + 1
        constraints%H_Group(1,constraints%nh(2),2) = i1
        constraints%H_Group(2,constraints%nh(2),2) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(2),2) = constraints%H_index(2,i1)

      else if (constraints%duplicate(i1) == 3) then

        constraints%nh(3) = constraints%nh(3) + 1
        constraints%H_Group(1,constraints%nh(3),3) = i1
        constraints%H_Group(2,constraints%nh(3),3) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(3),3) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(3),3) = constraints%H_index(3,i1)

      else if (constraints%duplicate(i1) == 4) then

        constraints%nh(4) = constraints%nh(4) + 1
        constraints%H_Group(1,constraints%nh(4),4) = i1
        constraints%H_Group(2,constraints%nh(4),4) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(4),4) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(4),4) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(4),4) = constraints%H_index(4,i1)

      else if (constraints%duplicate(i1) == 5) then

        constraints%nh(5) = constraints%nh(5) + 1
        constraints%H_Group(1,constraints%nh(5),5) = i1
        constraints%H_Group(2,constraints%nh(5),5) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(5),5) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(5),5) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(5),5) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(5),5) = constraints%H_index(5,i1)

      else if (constraints%duplicate(i1) == 6) then

        constraints%nh(6) = constraints%nh(6) + 1
        constraints%H_Group(1,constraints%nh(6),6) = i1
        constraints%H_Group(2,constraints%nh(6),6) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(6),6) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(6),6) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(6),6) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(6),6) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(6),6) = constraints%H_index(6,i1)

      else if (constraints%duplicate(i1) == 7) then

        constraints%nh(7) = constraints%nh(7) + 1
        constraints%H_Group(1,constraints%nh(7),7) = i1
        constraints%H_Group(2,constraints%nh(7),7) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(7),7) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(7),7) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(7),7) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(7),7) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(7),7) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(7),7) = constraints%H_index(7,i1)

      else if (constraints%duplicate(i1) == 8) then

        constraints%nh(8) = constraints%nh(8) + 1
        constraints%H_Group(1,constraints%nh(8),8) = i1
        constraints%H_Group(2,constraints%nh(8),8) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(8),8) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(8),8) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(8),8) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(8),8) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(8),8) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(8),8) = constraints%H_index(7,i1)
        constraints%H_Group(9,constraints%nh(8),8) = constraints%H_index(8,i1)

      end if

    end do


    return

  end subroutine setup_hbond_group_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_HBond_fep
  !> @brief        setup atom maps with H-bond connection groups for
  !                FEP calculations
  !! @authors      NK
  !! @param[in]    molecule    : molecule information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_HBond_fep(molecule, boundary, enefunc, constraints, &
                                 domain)

    ! formal arguments
    type(s_molecule),    target, intent(in)    :: molecule
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_enefunc),     target, intent(in)    :: enefunc
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(dp)                     :: x_shift, y_shift, z_shift
    real(dp)                     :: move(3), origin(3)
    integer                      :: i, j, icx, icy, icz, icel
    integer                      :: isolute, iwater, ih1, ih2, id
    integer                      :: icel_local
    integer                      :: ncel_local, ncel
    integer                      :: patm, psol, pwat, pnoh, phgl
    character(4)                 :: ci1
    logical                      :: mi1, cl1, tip4

    real(dp),            pointer :: bsize_x, bsize_y, bsize_z
    real(dp),            pointer :: csize_x, csize_y, csize_z
    integer,             pointer :: ncel_x, ncel_y, ncel_z
    integer,             pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: natom(:), nsolute(:), nwater(:)
    integer,             pointer :: solute_list(:,:), water_list(:,:,:)
    integer,             pointer :: No_HGr(:), HGr_local(:,:)
    integer,             pointer :: HGr_bond_list(:,:,:,:)

    ncel          = domain%num_cell_local + domain%num_cell_boundary
    ncel_local    = domain%num_cell_local

    call alloc_constraints(constraints, ConstraintsDomainBond_FEP, ncel, &
                           constraints%connect)

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z
    ncel_x        => boundary%num_cells_x
    ncel_y        => boundary%num_cells_y
    ncel_z        => boundary%num_cells_z
    csize_x       => boundary%cell_size_x
    csize_y       => boundary%cell_size_y
    csize_z       => boundary%cell_size_z

    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water
    solute_list   => domain%solute_list
    water_list    => domain%water_list

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list

    tip4          = constraints%tip4
    origin(1)     = boundary%origin_x
    origin(2)     = boundary%origin_y
    origin(3)     = boundary%origin_z


    ! solute atoms (not bonded to hydrogen) in each domain
    !
    do isolute = 1, enefunc%table%num_solute

      i   = enefunc%table%solute_list(isolute)
      ci1 = molecule%atom_name(i)

      mi1 = molecule%light_atom_mass(i)
      cl1 = molecule%light_atom_name(i)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1) 
      endif

!      if (constraints%duplicate(i) == 0 .and. ci1(1:1) /= 'H') then
      if (constraints%duplicate(i) == 0 .and. .not. cl1) then

        !coordinate shifted against the origin
        !
        x_shift = molecule%atom_coord(1,i) - boundary%origin_x
        y_shift = molecule%atom_coord(2,i) - boundary%origin_y
        z_shift = molecule%atom_coord(3,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        ! atoms inside the domain
        !
        if (cell_g2l(icel) /= 0) then

          ! local cell index
          !
          icel_local = cell_g2l(icel)

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          pnoh =  No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          call setup_pmelist_fep(molecule, i, &
                                  domain, icel_local, patm)

          natom  (icel_local) = patm
          nsolute(icel_local) = psol
          No_HGr (icel_local) = pnoh

        ! atoms in the boundary
        !
        else if (cell_g2b(icel) /= 0) then

          ! local cell index
          !
          icel_local = cell_g2b(icel) + ncel_local

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          pnoh =  No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          call setup_pmelist_fep(molecule, i, &
                                  domain, icel_local, patm)

          natom  (icel_local) = patm
          nsolute(icel_local) = psol
          No_HGr (icel_local) = pnoh

        end if 

      end if

    end do

    ! Hydrogen bonding group
    !
    do isolute = 1, constraints%connect

      do j = 1, constraints%nh(isolute)

        i = constraints%H_Group(1,j,isolute)
        x_shift = molecule%atom_coord(1,i) - boundary%origin_x
        y_shift = molecule%atom_coord(2,i) - boundary%origin_y
        z_shift = molecule%atom_coord(3,i) - boundary%origin_z

        move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1
        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        ! atoms inside the domain
        !
        if (cell_g2l(icel) /= 0) then
  
          icel_local = cell_g2l(icel)

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          phgl =  HGr_local(isolute,icel_local)
  
          ! hydrogen atom
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          call setup_pmelist_fep(molecule, i, &
                                  domain, icel_local, patm)

          HGr_local(isolute,icel_local) = phgl

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1,j,isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm

            call molecule_to_domain_fep(molecule, move, origin, i, &
                                    domain, icel_local, patm)

            call setup_pmelist_fep(molecule, i, &
                                    domain, icel_local, patm)

          end do

          natom  (icel_local) = patm
          nsolute(icel_local) = psol

        else if (cell_g2b(icel) /= 0) then

          icel_local = cell_g2b(icel) + ncel_local

          patm =  natom  (icel_local)
          psol =  nsolute(icel_local)
          phgl =  HGr_local(isolute,icel_local)

          ! hydrogen atoms
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, i, &
                                  domain, icel_local, patm)

          call setup_pmelist_fep(molecule, i, &
                                  domain, icel_local, patm)

          HGr_local(isolute,icel_local) = phgl

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1,j,isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm

            call molecule_to_domain_fep(molecule, move, origin, i, &
                                    domain, icel_local, patm)

            call setup_pmelist_fep(molecule, i, &
                                    domain, icel_local, patm)

          end do

          natom  (icel_local) = patm
          nsolute(icel_local) = psol

        end if

      end do

    end do

    ! water atoms in each domain
    !
    do iwater = 1, enefunc%table%num_water

      i   = enefunc%table%water_list(1,iwater)
      ih1 = enefunc%table%water_list(2,iwater)
      ih2 = enefunc%table%water_list(3,iwater)
      if (tip4) id = enefunc%table%water_list(4,iwater)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih2, &
                                domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, id, &
                                domain, icel_local, patm)

          call setup_pmelist_fep(molecule, id, &
                                domain, icel_local, patm)

        end if

        natom(icel_local)  = patm
        nwater(icel_local) = pwat

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih2, &
                                domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, id, &
                                domain, icel_local, patm)

          call setup_pmelist_fep(molecule, id, &
                                domain, icel_local, patm)

        end if

        natom(icel_local)  = patm
        nwater(icel_local) = pwat

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_HBond_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_table_fep
  !> @brief        setup atom maps with solute and water table
  !! @authors      NK
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_table_fep(molecule, boundary, enefunc, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3), origin(3)
    integer                   :: i, icx, icy, icz, icel
    integer                   :: isolute, iwater, ih1, ih2, id
    integer                   :: icel_local
    integer                   :: ncel_local, ncel
    integer                   :: patm, psol, pwat
    logical                   :: tip4

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: solute_list(:,:), water_list(:,:,:)

    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    ncel_x       => boundary%num_cells_x
    ncel_y       => boundary%num_cells_y
    ncel_z       => boundary%num_cells_z
    csize_x      => boundary%cell_size_x
    csize_y      => boundary%cell_size_y
    csize_z      => boundary%cell_size_z

    cell_g2l     => domain%cell_g2l
    cell_g2b     => domain%cell_g2b
    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water
    solute_list  => domain%solute_list
    water_list   => domain%water_list

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z

    tip4         = enefunc%table%tip4

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local

    ! solute atoms in each domain
    !
    do isolute = 1, enefunc%table%num_solute

      i = enefunc%table%solute_list(isolute)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)
      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom  (icel_local)
        psol =  nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1
        solute_list(psol,icel_local) = patm
        solute_list(patm,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        natom  (icel_local) = patm
        nsolute(icel_local) = psol

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom  (icel_local)
        psol =  nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1
        solute_list(psol,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        natom  (icel_local) = patm
        nsolute(icel_local) = psol

      end if

    end do

    ! water atoms in each domain
    !
    do iwater = 1, enefunc%table%num_water

      i   = enefunc%table%water_list(1,iwater)
      ih1 = enefunc%table%water_list(2,iwater)
      ih2 = enefunc%table%water_list(3,iwater)
      if (tip4) id = enefunc%table%water_list(4,iwater)

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih2, &
                                domain, icel_local, patm)

        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, id, &
                                domain, icel_local, patm)

          call setup_pmelist_fep(molecule, id, &
                                domain, icel_local, patm)

        end if

        natom (icel_local) = patm
        nwater(icel_local) = pwat

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm =  natom (icel_local)
        pwat =  nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih1, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih1, &
                                domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_fep(molecule, move, origin, ih2, &
                                domain, icel_local, patm)

        call setup_pmelist_fep(molecule, ih2, &
                                domain, icel_local, patm)

        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_fep(molecule, move, origin, id, &
                                domain, icel_local, patm)

          call setup_pmelist_fep(molecule, id, &
                                domain, icel_local, patm)

        end if

        natom (icel_local) = patm
        nwater(icel_local) = pwat

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_table_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_fep
  !> @brief        setup atom maps with whole atoms for Hybrid topology FEP calculation
  !! @authors      NK
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_fep(molecule, boundary, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    real(dp)                  :: x_shift, y_shift, z_shift
    real(dp)                  :: move(3), origin(3)
    integer                   :: i, icx, icy, icz, icel
    integer                   :: icel_local
    integer                   :: ncel_local, ncel, natom_all

    real(dp),         pointer :: bsize_x, bsize_y, bsize_z
    real(dp),         pointer :: csize_x, csize_y, csize_z
    integer,          pointer :: fep_factor(:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), nsolute(:)


    bsize_x        => boundary%box_size_x
    bsize_y        => boundary%box_size_y
    bsize_z        => boundary%box_size_z
    csize_x        => boundary%cell_size_x
    csize_y        => boundary%cell_size_y
    csize_z        => boundary%cell_size_z
    ncel_x         => boundary%num_cells_x
    ncel_y         => boundary%num_cells_y
    ncel_z         => boundary%num_cells_z

    cell_g2l          => domain%cell_g2l
    cell_g2b          => domain%cell_g2b
    natom             => domain%num_atom
    nsolute      => domain%num_solute

    origin(1)      = boundary%origin_x
    origin(2)      = boundary%origin_y
    origin(3)      = boundary%origin_z

    ncel           = domain%num_cell_local + domain%num_cell_boundary
    ncel_local     = domain%num_cell_local
    natom_all      = domain%num_atom_all


    do i = 1, natom_all

      !coordinate shifted against the origin
      !
      x_shift = molecule%atom_coord(1,i) - boundary%origin_x
      y_shift = molecule%atom_coord(2,i) - boundary%origin_y
      z_shift = molecule%atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_dp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_dp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_dp - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      !assign which cell
      !
      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y


      ! atoms inside the domain
      !
      if (cell_g2l(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2l(icel)

        natom(icel_local) = natom(icel_local) + 1

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, natom(icel_local))

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, natom(icel_local))

      ! atoms in a boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        natom(icel_local) = natom(icel_local) + 1

        call molecule_to_domain_fep(molecule, move, origin, i, &
                                domain, icel_local, natom(icel_local))

        call setup_pmelist_fep(molecule, i, &
                                domain, icel_local, natom(icel_local))

      end if


    end do

    nsolute(1:ncel) = natom(1:ncel)
    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain_fep
  !> @brief        copy molecule information to domain for FEP calculations
  !! @authors      NK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain_fep(molecule, move, origin, iatom, &
                                domain, icel, icel_atom)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    real(dp),                 intent(in)    :: move(3)
    real(dp),                 intent(in)    :: origin(3)
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    ! local variables
    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)

    real(dp),         pointer :: coord_local (:,:,:)
    real(dp),         pointer :: ref_local   (:,:,:)
    real(dp),         pointer :: vel_local   (:,:,:)
    real(wp),         pointer :: charge_local(:,:)
    real(dp),         pointer :: mass_local  (:,:)
    real(wp),         pointer :: trans       (:,:,:)
    integer,          pointer :: class_local (:,:)
    integer,          pointer :: id_l2g      (:,:)
    integer,          pointer :: id_g2l      (:,:)


    coord            => molecule%atom_coord
    vel              => molecule%atom_velocity
    charge           => molecule%charge
    mass             => molecule%mass
    atom_class       => molecule%atom_cls_no

    id_g2l       => domain%id_g2l
    id_l2g       => domain%id_l2g
    coord_local  => domain%coord
    ref_local    => domain%coord_ref
    vel_local    => domain%velocity
    charge_local => domain%charge
    mass_local   => domain%mass
    class_local  => domain%atom_cls_no
    trans        => domain%trans_vec


    id_l2g(icel_atom,icel) = iatom
    id_g2l(1,iatom) = icel
    id_g2l(2,iatom) = icel_atom


    coord_local(1:3,icel_atom,icel) = coord (1:3,iatom) - origin(1:3)
    ref_local  (1:3,icel_atom,icel) = coord (1:3,iatom) - origin(1:3)
    vel_local  (1:3,icel_atom,icel) = vel   (1:3,iatom)
    charge_local   (icel_atom,icel) = charge    (iatom)
    mass_local     (icel_atom,icel) = mass      (iatom)
    class_local    (icel_atom,icel) = atom_class(iatom)
    trans      (1:3,icel_atom,icel) = move(1:3) 
    domain%fepgrp  (icel_atom,icel) = int(molecule%fepgrp(iatom))

    return

  end subroutine molecule_to_domain_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord_fep
  !> @brief        check_atom_coordinate
  !! @authors      NK
  !! @param[in]    ene_info      : ENERGY section control parameters information
  !! @param[in]    boundary      : boundary condition information
  !! @param[in]    contact_check : flag for contact_check
  !! @param[inout] domain        : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord_fep(ene_info, boundary, contact_check, domain)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    type(s_boundary), target, intent(in)    :: boundary
    logical,                  intent(in)    :: contact_check
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    integer                   :: i, ix, iy, ij, j
    integer                   :: id, my_id,  omp_get_thread_num
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: cell_pairlist1(:,:)
    integer,          pointer :: id_l2g(:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(dp),         pointer :: mass(:,:), coord(:,:,:) 
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)

    ! FEP
    integer                   :: fg1, fg2

    mass           => domain%mass 
    coord          => domain%coord
    trans1         => domain%trans_vec
    trans2         => domain%translated
    cell_move      => domain%cell_move
    system_size    => domain%system_size
    natom          => domain%num_atom
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    cell_pairlist1 => domain%cell_pairlist1
    id_l2g         => domain%id_l2g

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
        end do
      end do
    else
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(1:3,ix,i) = coord(1:3,ix,i) 
        end do
      end do
    endif
    !$omp barrier

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        dir(1:3) = trans2(1:3,ix,i)
        do iy = ix + 1, natom(i)
          dr(1:3) = dir(1:3) - trans2(1:3,iy,i)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)

          fg1 = domain%fepgrp(ix,i)
          fg2 = domain%fepgrp(iy,i)
          if((fg1 == 1) .and. (fg2 == 2)) cycle
          if((fg1 == 2) .and. (fg2 == 1)) cycle
          if((fg1 == 1) .and. (fg2 == 4)) cycle
          if((fg1 == 4) .and. (fg2 == 1)) cycle
          if((fg1 == 2) .and. (fg2 == 3)) cycle
          if((fg1 == 3) .and. (fg2 == 2)) cycle
          if((fg1 == 3) .and. (fg2 == 4)) cycle
          if((fg1 == 4) .and. (fg2 == 3)) cycle
          if((fg1 == 5) .and. (fg2 == 3)) cycle
          if((fg1 == 3) .and. (fg2 == 5)) cycle
          if((fg1 == 5) .and. (fg2 == 4)) cycle
          if((fg1 == 4) .and. (fg2 == 5)) cycle
          if((fg1 == 3) .and. (fg2 == 3)) cycle
          if((fg1 == 4) .and. (fg2 == 4)) cycle
          if((fg1 == 1) .and. (fg2 == 3)) cycle
          if((fg1 == 3) .and. (fg2 == 1)) cycle
          if((fg1 == 2) .and. (fg2 == 4)) cycle
          if((fg1 == 4) .and. (fg2 == 2)) cycle

          if (rr < EPS) then
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
            call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
          end if

          if (rr < ene_info%minimum_contact .and. &
              abs(mass(ix,i)) > EPS .and. abs(mass(iy,i)) > EPS) then 
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, i), sqrt(rr)
            !$omp end critical
            if (rr < ene_info%err_minimum_contact .and.  &
               .not. contact_check) &
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
          endif
        end do
      end do
    end do

    if (boundary%type == BoundaryTypePBC) then
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)

        tr(1:3) = cell_move(1:3,j,i)*system_size(1:3)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i) 
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - coord(1:3,iy,j)+ tr(1:3)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)

            fg1 = domain%fepgrp(ix,i)
            fg2 = domain%fepgrp(iy,j)
            if((fg1 == 1) .and. (fg2 == 2)) cycle
            if((fg1 == 2) .and. (fg2 == 1)) cycle
            if((fg1 == 1) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 2)) cycle
            if((fg1 == 3) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 3)) cycle
            if((fg1 == 5) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 5)) cycle
            if((fg1 == 5) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 5)) cycle
            if((fg1 == 3) .and. (fg2 == 3)) cycle
            if((fg1 == 4) .and. (fg2 == 4)) cycle
            if((fg1 == 1) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 2)) cycle

            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < ene_info%minimum_contact .and. &
                abs(mass(ix,i)) > EPS .and. abs(mass(iy,j)) > EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and.  &
                 .not. contact_check) &
                call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
            endif
     
          end do
        end do
      end do
    else
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)
        do ix = 1, natom(i)
     
          dir(1:3) = trans2(1:3,ix,i)
     
          do iy = 1, natom(j)
            dr(1:3) = dir(1:3) - trans2(1:3,iy,j)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)

            fg1 = domain%fepgrp(ix,i)
            fg2 = domain%fepgrp(iy,j)
            if((fg1 == 1) .and. (fg2 == 2)) cycle
            if((fg1 == 2) .and. (fg2 == 1)) cycle
            if((fg1 == 1) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 2)) cycle
            if((fg1 == 3) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 3)) cycle
            if((fg1 == 5) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 5)) cycle
            if((fg1 == 5) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 5)) cycle
            if((fg1 == 3) .and. (fg2 == 3)) cycle
            if((fg1 == 4) .and. (fg2 == 4)) cycle
            if((fg1 == 1) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 2)) cycle

            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual)')
            end if

            if (rr < ene_info%minimum_contact .and. &
                abs(mass(ix,i)) > EPS .and. abs(mass(iy,j)) > EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
                'WARNING: too short distance:',id_l2g(ix,i), id_l2g(iy, j), sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and.  &
                 .not. contact_check) &
                call error_msg('Check_Atom_Coord> Some atoms have large clashes (see "Chapter: Trouble shooting" in the user manual). Please use nonb_limiter.')
            endif
     
          end do
        end do
      end do
    endif

    !$omp end parallel

    return

  end subroutine check_atom_coord_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pmelist_fep(molecule, iatom, &
                                domain, icel, icel_atom)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    ! local variables
    integer,             pointer :: natom_preserve(:)
    integer,             pointer :: natom_singleA(:)
    integer,             pointer :: natom_singleB(:)
    integer,             pointer :: natom_appear_gr(:)
    integer,             pointer :: natom_vanish_gr(:)
    integer,             pointer :: pmelist_preserve(:,:)
    integer,             pointer :: pmelist_appear_gr(:,:)
    integer,             pointer :: pmelist_vanish_gr(:,:)
    integer,             pointer :: id_singleA(:,:,:)
    integer,             pointer :: id_singleB(:,:,:)

    natom_preserve    => domain%num_atom_preserve
    natom_singleA     => domain%num_atom_singleA
    natom_singleB     => domain%num_atom_singleB
    natom_appear_gr   => domain%num_atom_appear_gr
    natom_vanish_gr   => domain%num_atom_vanish_gr
    pmelist_preserve  => domain%pmelist_preserve
    pmelist_appear_gr => domain%pmelist_appear_gr
    pmelist_vanish_gr => domain%pmelist_vanish_gr
    id_singleA        => domain%id_singleA
    id_singleB        => domain%id_singleB

    if(molecule%fep_topology == 2) then
      if (int(molecule%fepgrp(iatom)) == 5) then
        natom_preserve(icel)  = natom_preserve(icel)  + 1
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_preserve(natom_preserve(icel),icel)   = icel_atom
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 1) then
        natom_preserve(icel)  = natom_preserve(icel)  + 1
        natom_singleA(icel)   = natom_singleA(icel)   + 1
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_preserve(natom_preserve(icel),icel)   = icel_atom
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 2) then
        natom_preserve(icel)  = natom_preserve(icel)  + 1
        natom_singleB(icel)   = natom_singleB(icel)   + 1
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_preserve(natom_preserve(icel),icel)   = icel_atom
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 3) then
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 4) then
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
      end if
    else
      if (int(molecule%fepgrp(iatom)) == 5) then
        natom_preserve(icel)  = natom_preserve(icel)  + 1
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_preserve(natom_preserve(icel),icel)   = icel_atom
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 1) then
        natom_singleA(icel)   = natom_singleA(icel)   + 1
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
        id_singleA(natom_singleA(icel),icel,1) = icel_atom
        id_singleA(natom_singleA(icel),icel,2) = domain%id_l2g(icel_atom, icel)
        id_singleA(natom_singleA(icel),icel,3) = natom_singleA(icel)
      else if (int(molecule%fepgrp(iatom)) == 2) then
        natom_singleB(icel)   = natom_singleB(icel)   + 1
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
        id_singleB(natom_singleB(icel),icel,1) = icel_atom
        id_singleB(natom_singleB(icel),icel,2) = domain%id_l2g(icel_atom, icel)
        id_singleB(natom_singleB(icel),icel,3) = natom_singleB(icel)
      else if (int(molecule%fepgrp(iatom)) == 3) then
        natom_vanish_gr(icel) = natom_vanish_gr(icel) + 1
        pmelist_vanish_gr(natom_vanish_gr(icel),icel) = icel_atom
      else if (int(molecule%fepgrp(iatom)) == 4) then
        natom_appear_gr(icel) = natom_appear_gr(icel) + 1
        pmelist_appear_gr(natom_appear_gr(icel),icel) = icel_atom
      end if
    end if

    return

  end subroutine setup_pmelist_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sort_single_fep
  !> @brief        sort atoms in single region in FEP using heap sort
  !                in order to match the order of atoms in single region
  !! @authors      HO
  !! @param[inout] domain   : domain information
  !! @param[in]    enefunc  : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sort_single_fep(domain, enefunc)
    implicit none

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc

    ! local variables
    integer          :: i, j, k, l, icel
    integer          :: temp, temp2
    integer, pointer :: natom_singleA(:)
    integer, pointer :: natom_singleB(:)
    integer, pointer :: id_singleA(:,:,:)
    integer, pointer :: id_singleB(:,:,:)

    natom_singleA => domain%num_atom_singleA
    natom_singleB => domain%num_atom_singleB
    id_singleA    => domain%id_singleA
    id_singleB    => domain%id_singleB

    if (enefunc%fep_topology == 2) return

    ! Sort singleA
    !
    do icel = 1, domain%num_cell_local + domain%num_cell_boundary
      if (natom_singleA(icel) <= 1) cycle
      l = natom_singleA(icel)/2 + 1
      k = natom_singleA(icel)
      do while (.true.)
        if (l > 1) then
          l = l - 1
          temp  = id_singleA(l,icel,2)
          temp2 = id_singleA(l,icel,3)
        else
          temp  = id_singleA(k,icel,2)
          temp2 = id_singleA(k,icel,3)
          id_singleA(k,icel,2) = id_singleA(1,icel,2)
          id_singleA(k,icel,3) = id_singleA(1,icel,3)
          k = k - 1
          if (k == 1) then
            id_singleA(1,icel,2) = temp
            id_singleA(1,icel,3) = temp2
            exit
          end if
        end if
        i = l
        j = 2 * l
        do while (j <= k)
          if (j < k) then
            if (id_singleA(j,icel,2) < id_singleA(j+1,icel,2)) then
              j = j + 1
            end if
          end if
          if (temp < id_singleA(j,icel,2)) then
            id_singleA(i,icel,2) = id_singleA(j,icel,2)
            id_singleA(i,icel,3) = id_singleA(j,icel,3)
            i = j
            j = 2 * j
          else
            exit
          end if
        end do
        id_singleA(i,icel,2) = temp
        id_singleA(i,icel,3) = temp2
      end do
    end do

    ! Sort singleB
    !
    do icel = 1, domain%num_cell_local + domain%num_cell_boundary
      if (natom_singleB(icel) <= 1) cycle
      l = natom_singleB(icel)/2 + 1
      k = natom_singleB(icel)
      do while (.true.)
        if (l > 1) then
          l = l - 1
          temp  = id_singleB(l,icel,2)
          temp2 = id_singleB(l,icel,3)
        else
          temp  = id_singleB(k,icel,2)
          temp2 = id_singleB(k,icel,3)
          id_singleB(k,icel,2) = id_singleB(1,icel,2)
          id_singleB(k,icel,3) = id_singleB(1,icel,3)
          k = k - 1
          if (k == 1) then
            id_singleB(1,icel,2) = temp
            id_singleB(1,icel,3) = temp2
            exit
          end if
        end if
        i = l
        j = 2 * l
        do while (j <= k)
          if (j < k) then
            if (id_singleB(j,icel,2) < id_singleB(j+1,icel,2)) then
              j = j + 1
            end if
          end if
          if (temp < id_singleB(j,icel,2)) then
            id_singleB(i,icel,2) = id_singleB(j,icel,2)
            id_singleB(i,icel,3) = id_singleB(j,icel,3)
            i = j
            j = 2 * j
          else
            exit
          end if
        end do
        id_singleB(i,icel,2) = temp
        id_singleB(i,icel,3) = temp2
      end do
    end do

    return

  end subroutine sort_single_fep

end module sp_domain_mod
