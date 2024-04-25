!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_update_domain_mod
!> @brief   library subroutine used for integrator
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_update_domain_mod

  use sp_domain_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_communicate_mod
  use sp_migration_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: domain_interaction_update_md
  public  :: domain_interaction_update_min
  public  :: do_migration_nobc
  !FEP
  public  :: domain_interaction_update_md_fep
  public  :: domain_interaction_update_min_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_md
  !> @brief        update interactions with MD
  !! @authors      JJ
  !! @param[in]    istep       : step number
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] constraints : constraints information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_md(istep, dynamics, domain, enefunc, &
                                          pairlist, boundary, constraints, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0

    if (constraints%tip4 .or. enefunc%table%tip4) then
      water_atom = 4
    else
      water_atom = 3
    end if

    if (mod(istep,dynamics%nbupdate_period) == 0) then

#ifdef DEBUG
      if (main_rank) then
        write(MsgOut,'(a,3f15.8)') &
             "Debuging > bsize_[x,y,z]", &
             boundary%box_size_x,boundary%box_size_y,boundary%box_size_z
        write(MsgOut,'(a,2i8)')    &
             "Debuging > ncell, nb",ncell,nb
      end if
#endif
    end if

    if (mod(istep,dynamics%nbupdate_period) == 0) then

      if (boundary%type == BoundaryTypeNOBC) then

        call update_nobc_boundary(domain, boundary)
        if (istep == dynamics%nbupdate_period) then
          call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
          call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)
        end if
        call do_migration_nobc(domain, enefunc, boundary, comm, pairlist)

      else

        pairlist%univ_update = 1

        if (constraints%rigid_bond) then

          if (istep == dynamics%nbupdate_period) then

            call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
            call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)

            call alloc_constraints(constraints, ConstraintsHGroupMove, &
                                   constraints%connect, ncell+nb)

            call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)

          end if

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerUpdate, TimerOn)
          call update_outgoing_atom (boundary, constraints, domain)
          call update_outgoing_HGr  (boundary, constraints, domain)
          call update_outgoing_water(water_atom, boundary, domain)
          call timer(TimerUpdate, TimerOff)
  
          call communicate_constraints(domain, comm, constraints)
  
          call timer(TimerUpdate, TimerOn)
          call update_incoming_atom (constraints, domain)
          call update_incoming_HGr  (constraints, domain)
          call update_incoming_water(water_atom, domain)
  
          call update_cell_size_constraints    (domain, comm, constraints)
          call update_cell_boundary_constraints(domain, comm, boundary, &
                                                constraints)
  
          call update_enefunc(.false., domain, comm, enefunc, constraints)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerIntegrator, TimerOff)
          call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
          if (enefunc%pairlist_check) then
            call update_pairlist_pbc_check(enefunc, domain, pairlist)
          else
            call update_pairlist_pbc(enefunc, domain, pairlist)
          endif
#else
          call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
          call timer(TimerPairList, TimerOff)

        else

          if (istep == dynamics%nbupdate_period) then

            call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
            call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)
            call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)
  
          end if

          call timer(TimerUpdate, TimerOn)
          call update_outgoing_solute (boundary, domain)
          call update_outgoing_water  (water_atom, boundary, domain)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerComm3, TimerOn)
          call communicate_solutewater(enefunc, domain, comm)
          call timer(TimerComm3, TimerOff)
  
          call timer(TimerUpdate, TimerOn)
          call update_incoming_solute (domain)
          call update_incoming_water  (water_atom, domain)
  
          call update_cell_size_solute   (domain, comm)
          call update_cell_boundary_table(enefunc, domain, comm, boundary)
  
          call update_enefunc(.false., domain, comm, enefunc)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerIntegrator, TimerOff)
          call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
          if (enefunc%pairlist_check) then
            call update_pairlist_pbc_check(enefunc, domain, pairlist)
          else
            call update_pairlist_pbc(enefunc, domain, pairlist)
          endif
#else
          call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
          call timer(TimerPairList, TimerOff)

        end if

      end if
    end if

    return

  end subroutine domain_interaction_update_md
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_min
  !> @brief        update interactions with minimization
  !! @authors      JJ
  !! @param[in]    istep       : step number
  !! @param[in]    minimize    : minimize information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_min(istep, minimize, domain, enefunc, &
                                           pairlist, boundary, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_minimize),        intent(in)    :: minimize
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom

    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0
    if (enefunc%table%tip4) then
      water_atom = 4
    else
      water_atom = 3
    end if

    if (mod(istep,minimize%nbupdate_period) == 0) then

      pairlist%univ_update = 1
      call timer(TimerIntegrator, TimerOn)

      if (istep == minimize%nbupdate_period) then

        call alloc_domain(domain, DomainPtlMove,   ncell+nb, ncell, 1)
        call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)
        call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)

      end if

      call update_outgoing_solute (boundary, domain)
      call update_outgoing_water  (water_atom, boundary, domain)

      call communicate_solutewater(enefunc, domain, comm)

      call update_incoming_solute (domain)
      call update_incoming_water  (water_atom, domain)

      call update_cell_size_solute   (domain, comm)
      call update_cell_boundary_table(enefunc, domain, comm, boundary)

      call update_enefunc(.false., domain, comm, enefunc)

      call timer(TimerIntegrator, TimerOff)
      call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc(enefunc, domain, pairlist)
      endif
#else
      call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine domain_interaction_update_min

  subroutine do_migration_nobc(domain, enefunc, boundary, comm, pairlist)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_boundary),        intent(inout) :: boundary
    type(s_comm),            intent(inout) :: comm
    type(s_pairlist),        intent(inout) :: pairlist

    call timer(TimerUpdate, TimerOn)
    call update_outgoing_ptl(boundary, domain)
    call timer(TimerUpdate, TimerOff)
    call timer(TimerComm3, TimerOn)
    call communicate_ptl    (domain, comm)
    call timer(TimerComm3, TimerOff)
    call timer(TimerUpdate, TimerOn)
    call update_incoming_ptl(domain)

    call update_cell_size    (domain, comm)
    call update_cell_boundary(domain, comm, boundary)

    call update_enefunc_contact(domain, comm, enefunc)
    call timer(TimerUpdate, TimerOff)

    call timer(TimerIntegrator, TimerOff)
    call timer(TimerPairList, TimerOn)
    call update_pairlist_nobc(enefunc, domain, pairlist)
    call timer(TimerPairList, TimerOff)

    return

  end subroutine do_migration_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_md_fep
  !> @brief        update interactions with MD for FEP calculations
  !! @authors      NK, HO
  !! @param[in]    istep       : step number
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] constraints : constraints information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_md_fep(istep, dynamics, domain, &
                                          enefunc, pairlist, boundary, &
                                          constraints, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0

    if (constraints%tip4 .or. enefunc%table%tip4) then
      water_atom = 4
    else
      water_atom = 3
    end if

    if (mod(istep,dynamics%nbupdate_period) == 0) then

#ifdef DEBUG
      if (main_rank) then
        write(MsgOut,'(a,3f15.8)') &
             "Debuging > bsize_[x,y,z]", &
             boundary%box_size_x,boundary%box_size_y,boundary%box_size_z
        write(MsgOut,'(a,2i8)')    &
             "Debuging > ncell, nb",ncell,nb
      end if
#endif
    end if

    if (mod(istep,dynamics%nbupdate_period) == 0) then

      if (boundary%type == BoundaryTypeNOBC) then

        call update_nobc_boundary(domain, boundary)
        if (istep == dynamics%nbupdate_period) then
          call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
          call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)
        end if
        call do_migration_nobc(domain, enefunc, boundary, comm, pairlist)

      else

        pairlist%univ_update = 1

        if (constraints%rigid_bond) then

          if (istep == dynamics%nbupdate_period) then

            call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
            call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)

            call alloc_constraints(constraints, ConstraintsHGroupMove, &
                                   constraints%connect, ncell+nb)

            call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)

          end if

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerUpdate, TimerOn)
          call update_outgoing_atom_fep (boundary, constraints, domain)
          call update_outgoing_HGr_fep  (boundary, constraints, domain)
          call update_outgoing_water_fep(water_atom, boundary, domain)
          call timer(TimerUpdate, TimerOff)
  
          call communicate_constraints_fep(domain, comm, constraints)
  
          call timer(TimerUpdate, TimerOn)
          call update_incoming_atom_fep (constraints, domain)
          call update_incoming_HGr_fep  (constraints, domain)
          call update_incoming_water_fep(water_atom, domain)
  
          call update_cell_size_constraints    (domain, comm, constraints)
          call update_cell_boundary_constraints_fep(domain, comm, boundary, &
                                                constraints)
  
          call update_pmelist_fep(domain, enefunc)
          call sort_single_fep(domain, enefunc)

          call update_enefunc(.false., domain, comm, enefunc, constraints)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerIntegrator, TimerOff)
          call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
          if (enefunc%pairlist_check) then
            call update_pairlist_pbc_check_fep(enefunc, domain, pairlist)
          else
            call update_pairlist_pbc_fep(enefunc, domain, pairlist)
          endif
#else
          call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
          call timer(TimerPairList, TimerOff)

        else

          if (istep == dynamics%nbupdate_period) then

            call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
            call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)
            call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)
  
          end if

          call timer(TimerUpdate, TimerOn)
          call update_outgoing_solute_fep (boundary, domain)
          call update_outgoing_water_fep  (water_atom, boundary, domain)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerComm3, TimerOn)
          call communicate_solutewater_fep(enefunc, domain, comm)
          call timer(TimerComm3, TimerOff)
  
          call timer(TimerUpdate, TimerOn)
          call update_incoming_solute_fep (domain)
          call update_incoming_water_fep  (water_atom, domain)
  
          call update_cell_size_solute   (domain, comm)
          call update_cell_boundary_table_fep(enefunc, domain, comm, boundary)
  
          call update_pmelist_fep(domain, enefunc)
          call sort_single_fep(domain, enefunc)

          call update_enefunc(.false., domain, comm, enefunc)
          call timer(TimerUpdate, TimerOff)
  
          call timer(TimerIntegrator, TimerOff)
          call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
          if (enefunc%pairlist_check) then
            call update_pairlist_pbc_check_fep(enefunc, domain, pairlist)
          else
            call update_pairlist_pbc_fep(enefunc, domain, pairlist)
          endif
#else
          call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
          call timer(TimerPairList, TimerOff)

        end if

      end if
    end if

    return

  end subroutine domain_interaction_update_md_fep
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update_min_fep
  !> @brief        update interactions with minimization for FEP
  !! @authors      HO, NK
  !! @param[in]    istep       : step number
  !! @param[in]    minimize    : minimize information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update_min_fep(istep, minimize, domain, enefunc, &
                                           pairlist, boundary, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    type(s_minimize),        intent(in)    :: minimize
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom

    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0
    if (enefunc%table%tip4) then
      water_atom = 4
    else
      water_atom = 3
    end if

    if (mod(istep,minimize%nbupdate_period) == 0) then

      pairlist%univ_update = 1
      call timer(TimerIntegrator, TimerOn)

      if (istep == minimize%nbupdate_period) then

        call alloc_domain(domain, DomainPtlMove,   ncell+nb, ncell, 1)
        call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, water_atom)
        call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)

      end if

      call update_outgoing_solute_fep (boundary, domain)
      call update_outgoing_water_fep  (water_atom, boundary, domain)

      call communicate_solutewater_fep(enefunc, domain, comm)

      call update_incoming_solute_fep (domain)
      call update_incoming_water_fep  (water_atom, domain)

      call update_cell_size_solute   (domain, comm)
      call update_cell_boundary_table_fep(enefunc, domain, comm, boundary)

      call update_pmelist_fep(domain, enefunc)
      call sort_single_fep(domain, enefunc)

      call update_enefunc(.false., domain, comm, enefunc)

      call timer(TimerIntegrator, TimerOff)
      call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check_fep(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc_fep(enefunc, domain, pairlist)
      endif
#else
      call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine domain_interaction_update_min_fep

end module sp_update_domain_mod
