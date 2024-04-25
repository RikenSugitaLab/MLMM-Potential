!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_setup_spdyn_peer_mod
!> @brief   peer module of lib/domain_setup_md_mod
!! @authors Norio Takase (NT), Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_setup_spdyn_peer_mod

  use pr_select_atoms_peer_mod
  use pr_control_mod
  use pr_domain_index_mod
  use pr_huge_molecule_mod
  use sp_dynamics_mod
  use sp_dynvars_mod
  use sp_domain_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_enefunc_restraints_mod
  use sp_energy_mod
  use sp_boundary_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use dihedral_libs_mod
  use select_mod
  use select_atoms_str_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: setup_md_peer_0
  private ::   setup_boundary_peer_0
  private ::   setup_domain_peer_0
  private ::     setup_solute_and_water_peer_0
  private ::     setup_hbond_group_peer_0
  private ::   setup_restraints_peer_0
  private ::   define_enefunc_peer_0
  private ::     define_enefunc_charmm_peer_0
  private ::       setup_enefunc_nonb_C_peer_0
  private ::     define_enefunc_amber_peer_0
  private ::       setup_enefunc_nonb_A_peer_0
  private ::     define_enefunc_gromacs_peer_0
  private ::       setup_enefunc_nonb_G_peer_0
  private ::       setup_enefunc_table_peer_0
  private ::       setup_enefunc_restraints_peer_0
  private ::         setup_enefunc_rest_group_peer_0
  private ::         setup_enefunc_rest_func_peer_0
  private ::       setup_enefunc_gro_restraints_peer_0
  private ::   setup_dynamics_peer_0
  private ::   setup_constraints_peer_0
  private ::     setup_fast_water_tip4_peer_0
  private ::     setup_fast_water_peer_0
  private ::     setup_rigid_bond_peer_0

  public  :: setup_md_peer
  private ::   setup_domain_peer
  private ::     setup_atom_by_HBond_peer
  private ::     setup_atom_by_table_peer
  private ::     setup_atom_peer
  private ::       molecule_to_domain_peer
  private ::   define_enefunc_peer
  private ::     define_enefunc_charmm_peer
  private ::       setup_enefunc_bond_C_peer
  private ::       setup_enefunc_angl_C_peer
  private ::       setup_enefunc_dihe_C_peer
  private ::       setup_enefunc_impr_C_peer
  private ::       setup_enefunc_cmap_C_peer
  private ::       setup_enefunc_nonb_C_peer
  private ::       setup_enefunc_bond_constraint_C_peer
  private ::       setup_enefunc_angl_constraint_C_peer
  private ::     define_enefunc_amber_peer
  private ::       setup_enefunc_bond_A_peer
  private ::       setup_enefunc_angl_A_peer
  private ::       setup_enefunc_dihe_A_peer
  private ::       setup_enefunc_impr_A_peer
  private ::       setup_enefunc_nonb_A_peer
  private ::       setup_enefunc_bond_constraint_A_peer
  private ::       setup_enefunc_angl_constraint_A_peer
  private ::     define_enefunc_gromacs_peer
  private ::       setup_enefunc_bond_G_peer
  private ::       setup_enefunc_angl_G_peer
  private ::       setup_enefunc_dihe_G_peer
  private ::       setup_enefunc_rb_dihe_G_peer
  private ::       setup_enefunc_impr_G_peer
  private ::       setup_enefunc_bond_constraint_G_peer
  private ::       setup_enefunc_angl_constraint_G_peer
  private ::       setup_enefunc_nonb_G_peer
  private ::       setup_enefunc_restraints_peer
  private ::         setup_enefunc_rest_domain_peer
  private ::     alloc_id_g2l
  
contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_md_peer_0
  !> @brief        peer routine of sp_setup_spdyn_mod::setup_md()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_md_peer_0(ctrl_data,    &
                             par,          &
                             prmtop,       &
                             grotop,       &
                             boundary_rst, &
                             boundary,     &
                             enefunc,      &
                             dynvars,      &
                             dynamics,     &
                             constraints,  &
                             restraints)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_boundary),        intent(in)    :: boundary_rst
    type(s_boundary),        intent(inout) :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(inout) :: restraints


    ! set parameters for boundary condition
    !
    call setup_boundary_peer_0(ctrl_data%bou_info,              &
                               ctrl_data%ene_info%table,        &
                               ctrl_data%ene_info%pairlistdist, &
                               ctrl_data%ene_info%water_model,  &
                               ctrl_data%ens_info%ensemble,     &
                               ctrl_data%cons_info%rigid_bond,  &
                               ctrl_data%ene_info%dsize_cg,     &
                               ctrl_data%ene_info%dmin_size_cg, &
                               boundary_rst,                    &
                               boundary)

    ! set parameters for domain 
    !
    call setup_domain_peer_0(ctrl_data%ene_info,  &
                             ctrl_data%cons_info, &
                             boundary, enefunc, constraints)


    ! set parameters for communication
    !
    ! ..skip


    ! set parameters for restraints
    !
    call setup_restraints_peer_0(ctrl_data%res_info, &
                                 ctrl_data%sel_info, &
                                 restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc_peer_0(ctrl_data%ene_info, &
                               par, prmtop, grotop, restraints, enefunc)


    ! set parameters for pairlist
    !
    ! ..
    ! ..skip


    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars)


    ! set parameters for dynamics
    !
    call setup_dynamics_peer_0(ctrl_data%dyn_info, ctrl_data%bou_info, &
                               ctrl_data%res_info, dynamics)


    ! set parameters for ensemble
    !
    ! ..skip


    ! set parameters for constraint
    !
    call setup_constraints_peer_0(ctrl_data%cons_info, &
                                     par, prmtop, grotop, enefunc, constraints)


    return

  end subroutine setup_md_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_peer_0
  !> @brief        peer routine of sp_boundary_mod :: setup_boundary()
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_peer_0(bound_info, table, pairlistdist, watermodel,&
                            ensemble, rigid_bond, dsize_cg, dmin_size_cg,  &
                            boundary_rst, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: watermodel
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    logical,                 intent(in)    :: dsize_cg
    real(wp),                intent(in)    :: dmin_size_cg
    type(s_boundary),        intent(in)    :: boundary_rst
    type(s_boundary),        intent(inout) :: boundary


    call init_boundary(boundary)

    if (.not.((bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 2 .and. &
               bound_info%domain_z == 1))) then

      if (bound_info%domain_x == 1 .or. &
          bound_info%domain_y == 1 .or. &
          bound_info%domain_z == 1 ) then
        call error_msg('Setup_Boundary> other than (2,1,1)/(2,2,1), '//&
                       'domain[x,y,z] should be larger than 1')
      end if

    end if

    boundary%type           = bound_info%type
    boundary%origin_x       = bound_info%origin_x
    boundary%origin_y       = bound_info%origin_y
    boundary%origin_z       = bound_info%origin_z
    boundary%box_size_x     = bound_info%pbc_info%box_size_x
    boundary%box_size_y     = bound_info%pbc_info%box_size_y
    boundary%box_size_z     = bound_info%pbc_info%box_size_z
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    
    if (boundary_rst%box_size_x /= 0.0_wp) then

      boundary%box_size_x     = boundary_rst%box_size_x
      boundary%box_size_y     = boundary_rst%box_size_y
      boundary%box_size_z     = boundary_rst%box_size_z
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

    end if

    call setup_processor_number(bound_info, &
                                table, pairlistdist, watermodel, &
                                ensemble, rigid_bond, boundary)

    call setup_boundary_cell   (table, pairlistdist, watermodel, &
                                ensemble, rigid_bond, dsize_cg, dmin_size_cg,  &
                                boundary)

    return

  end subroutine setup_boundary_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_peer_0
  !> @brief        peer routine of sp_domain_mod::setup_domain()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_peer_0(ene_info, con_info, &
                                 boundary, enefunc, constraints)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    real(wp)                 :: mass_upper_bound

    ! initialize structure informations
    !
    call init_enefunc(enefunc)
    call init_constraints(constraints)

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type
    constraints%hydrogen_mass_upper_bound = &
                                con_info%hydrogen_mass_upper_bound
    mass_upper_bound = con_info%hydrogen_mass_upper_bound
    constraints%tip4 = .false.
    enefunc%table%tip4 = .false.

    if (constraints%rigid_bond) then

      call setup_solute_and_water_peer_0(mass_upper_bound,        &
                                         constraints%water_model, &
                                         constraints%tip4, enefunc)

      call setup_hbond_group_peer_0     (enefunc, constraints)

    else 

      call setup_solute_and_water_peer_0(mass_upper_bound,          &
                                         enefunc%table%water_model, &
                                         enefunc%table%tip4, enefunc)

    end if

    return

  end subroutine setup_domain_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_and_water_peer_0
  !> @brief        peer routine of sp_domain_mod::setup_solute_and_water()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_and_water_peer_0(mass_upper_bound, water_model, &
                                           tip4, enefunc)

    ! formal arguments
    real(wp),                intent(in)    :: mass_upper_bound
    character(*),            intent(in)    :: water_model
    logical,                 intent(inout) :: tip4
    type(s_enefunc),         intent(inout) :: enefunc
    
    ! local variables
    integer                  :: i, k, i1, i2, i3, io, ih, id(4)
    integer                  :: nwat, nsol, mole_no, tmp_mole_no
    integer                  :: water_atom
    character(4)             :: a1, a2, a3
    character(6)             :: res
    real(wp)                 :: mass(4), mass_max, mass_min


    ! s_enefunc%table
    !   num_water
    !   num_solute
    !   water_list
    !   solute_list
    !   solute_list_inv
    !

    nwat = 0
    nsol = 0

    do i = 1, hm_num_atoms
      res     = hm_residue_name(i)
      a1      = hm_atom_name(i)
      mass(1) = hm_mass(i)
      if ((res(1:3) == 'TIP' .or. res(1:3) == 'WAT' .or.   &
           res(1:3) == 'SOL' .or. res(1:3) == 'SPC') .and. & 
           a1(1:1) /= 'H' .and. a1(1:1) /= 'h' .and.       &
           mass(1) > mass_upper_bound) then
        nwat = nwat + 1
      else if (res(1:3) /= 'TIP' .and. res(1:3) /= 'WAT' .and. &
               res(1:3) /= 'SOL' .and. res(1:3) /= 'SPC') then
        nsol = nsol + 1
      end if
    end do

    if ((hm_num_atoms-nsol)/4 == nwat .and. nwat > 0) then
      tip4 = .true.
      water_atom = 4
    else if ((hm_num_atoms-nsol)/3 == nwat) then
      tip4 = .false.
      water_atom = 3
    else
      call error_msg('Setup_Solute_And_Water_Peer_0> # of water is incorrect.')
    end if

    enefunc%table%num_water  = nwat
    enefunc%table%num_solute = nsol

    if (main_rank) &
      write(MsgOut,'(a,i9)') &
         'Setup_Solute_And_Water_Peer_0> # of water molecule    : ', nwat

    call hm_alloc_water_list(nwat, water_atom)
    call hm_alloc_solute_list(nsol)
    call hm_alloc_solute_list_inv(hm_num_atoms)

    do i = 1, hm_num_atoms
      call hm_set_solute_list_inv(i, 0)
    end do

    i = 1
    nwat = 0
    nsol = 0
    do while(.true.)

      if (i > hm_num_atoms) &
        exit

      res = hm_residue_name(i)

      if (res(1:3) == 'TIP' .or. res(1:3) == 'WAT' .or. &
          res(1:3) == 'SOL' .or. res(1:3) == 'SPC') then

        do k = 1, water_atom
          mass(k) = hm_mass(i-1+k)
        end do

        mass_max = -1000.0_wp
        mass_min =  1000.0_wp
        do k = 1, water_atom
          mass_max = max(mass_max, mass(k))
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

        nwat = nwat + 1
        do k = 1, water_atom
          call hm_set_water_list(k,nwat,id(k))
        end do
        i = i + water_atom

      else

        nsol = nsol + 1
        call hm_set_solute_list    (nsol, i)
        call hm_set_solute_list_inv(i, nsol)
        i = i + 1

      end if

    end do

    if (nwat /= hm_num_water_list) &
      call error_msg('Setup_Water> number of water is incorrect')

    if (main_rank) &
      write(MsgOut,'(a,i9)') &
         'Setup_Solute_And_Water_Peer_0> # of solute            : ', nsol


    if (main_rank) &
      write(MsgOut,'(a)') ' '

    ! set water oxygen and hydrogen
    !
    io = hm_water_list(1,1)
    ih = hm_water_list(2,1)

    enefunc%table%atom_cls_no_O = hm_atom_cls_no(io)
    enefunc%table%atom_cls_no_H = hm_atom_cls_no(ih)
    enefunc%table%charge_O      = hm_charge(io)
    enefunc%table%charge_H      = hm_charge(ih)
    enefunc%table%mass_O        = hm_mass(io)
    enefunc%table%mass_H        = hm_mass(ih)

    if (tip4) then
      io = hm_water_list(4,1)
      enefunc%table%atom_cls_no_D = hm_atom_cls_no(io)
      enefunc%table%charge_D      = hm_charge(io)
      enefunc%table%mass_D        = hm_mass(io)
    end if

    return

  end subroutine setup_solute_and_water_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hbond_group_peer_0
  !> @brief        peer routine of sp_domain_mod::setup_hbond_group()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hbond_group_peer_0(enefunc, constraints)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, l, i1, i2, idu
    integer                  :: natom, nbond, nsolute, nhmax
    character(6)             :: ci1, ci2
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2
    real(wp)                 :: mass


    nsolute = hm_num_solute_list

    nbond = hm_num_bonds
    natom = hm_num_atoms

    call hm_alloc_duplicate(natom)
    call hm_alloc_H_index  (natom)

    do i = 1, natom
      call hm_set_duplicate(i, 0)
      do j = 1, 8
        call hm_set_H_index(j, i, 0)
      end do
    end do

    do i = 1, nbond

      i1 = hm_bond_list(1,i)
      i2 = hm_bond_list(2,i)

      if (hm_solute_list_inv(i1) /= 0 .and. &
          hm_solute_list_inv(i2) /= 0) then

        mass = hm_mass(i1) 
        if (mass > EPS .and. mass < constraints%hydrogen_mass_upper_bound) then
          mi1 = .true.
        else
          mi1 = .false.
        end if

        mass = hm_mass(i2) 
        if (mass > EPS .and. mass < constraints%hydrogen_mass_upper_bound) then
          mi2 = .true.
        else
          mi2 = .false.
        end if

        cl1 = .false. 
        ci1 = hm_atom_name(i1)
        if (ci1(1:1) == 'H' .or. ci1(1:1) == 'h' .or. &
            ci1(1:1) == 'D' .or. ci1(1:1) == 'd') then
          cl1 = .true.
        end if
        ci1 = hm_atom_cls_name(i1)
        if (ci1(1:1) == 'H' .or. ci1(1:1) == 'h' .or. &
            ci1(1:1) == 'D' .or. ci1(1:1) == 'd') then
          cl1 = .true.
        end if

        cl2 = .false. 
        ci2 = hm_atom_name(i2)
        if (ci2(1:1) == 'H' .or. ci2(1:1) == 'h' .or. &
            ci2(1:1) == 'D' .or. ci2(1:1) == 'd') then
          cl2 = .true.
        end if
        ci2 = hm_atom_cls_name(i2)
        if (ci2(1:1) == 'H' .or. ci2(1:1) == 'h' .or. &
            ci2(1:1) == 'D' .or. ci2(1:1) == 'd') then
          cl2 = .true.
        end if

        if (constraints%hydrogen_type == ConstraintAtomMass) then
          cl1 = mi1
          cl2 = mi2
        else if (constraints%hydrogen_type == ConstraintAtomBoth) then
          cl1 = (cl1 .or. mi1)
          cl2 = (cl2 .or. mi2)
        end if

        if (cl1 .or. cl2) then

          if (cl1) then
            idu = hm_duplicate(i2) + 1
            call hm_set_duplicate(i2, idu)
            call hm_set_H_index  (idu, i2, i1)
          else
            idu = hm_duplicate(i1) + 1
            call hm_set_duplicate(i1, idu)
            call hm_set_H_index  (idu, i1, i2)
          end if

        end if

      end if

    end do

    ! count XHn group for each number n
    !
    constraints%nh(1:8) = 0
    do i = 1, nsolute
      i1 = hm_solute_list(i)
      if (hm_duplicate(i1) == 1) then
        constraints%nh(1) = constraints%nh(1) + 1
      else if (hm_duplicate(i1) == 2) then
        constraints%nh(2) = constraints%nh(2) + 1
      else if (hm_duplicate(i1) == 3) then
        constraints%nh(3) = constraints%nh(3) + 1
      else if (hm_duplicate(i1) == 4) then
        constraints%nh(4) = constraints%nh(4) + 1
      else if (hm_duplicate(i1) == 5) then
        constraints%nh(5) = constraints%nh(5) + 1
      else if (hm_duplicate(i1) == 6) then
        constraints%nh(6) = constraints%nh(6) + 1
      else if (hm_duplicate(i1) == 7) then
        constraints%nh(7) = constraints%nh(7) + 1
      else if (hm_duplicate(i1) == 8) then
        constraints%nh(8) = constraints%nh(8) + 1
      else if (hm_duplicate(i1) >= 8) then
        call error_msg("Bond(>8) for one atom is not considered")
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
    do i = 1, nsolute

      i1 = hm_solute_list(i)

      if (hm_duplicate(i1) == 1) then
    
        constraints%nh(1) = constraints%nh(1) + 1
        constraints%H_Group(1, constraints%nh(1), 1)= i1
        constraints%H_Group(2, constraints%nh(1), 1)= hm_H_index(1, i1)

      else if (hm_duplicate(i1) == 2) then

        constraints%nh(2) = constraints%nh(2) + 1
        constraints%H_Group(1, constraints%nh(2), 2)= i1       
        constraints%H_Group(2, constraints%nh(2), 2)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(2), 2)= hm_H_index(2, i1)

      else if (hm_duplicate(i1) == 3) then

        constraints%nh(3) = constraints%nh(3) + 1
        constraints%H_Group(1, constraints%nh(3), 3)= i1       
        constraints%H_Group(2, constraints%nh(3), 3)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(3), 3)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(3), 3)= hm_H_index(3, i1)

      
      else if (hm_duplicate(i1) == 4) then

        constraints%nh(4) = constraints%nh(4) + 1
        constraints%H_Group(1, constraints%nh(4), 4)= i1       
        constraints%H_Group(2, constraints%nh(4), 4)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(4), 4)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(4), 4)= hm_H_index(3, i1)
        constraints%H_Group(5, constraints%nh(4), 4)= hm_H_index(4, i1)


      else if (hm_duplicate(i1) == 5) then

        constraints%nh(5) = constraints%nh(5) + 1
        constraints%H_Group(1, constraints%nh(5), 5)= i1       
        constraints%H_Group(2, constraints%nh(5), 5)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(5), 5)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(5), 5)= hm_H_index(3, i1)
        constraints%H_Group(5, constraints%nh(5), 5)= hm_H_index(4, i1)
        constraints%H_Group(6, constraints%nh(5), 5)= hm_H_index(5, i1)


      else if (hm_duplicate(i1) == 6) then

        constraints%nh(6) = constraints%nh(6) + 1
        constraints%H_Group(1, constraints%nh(6), 6)= i1       
        constraints%H_Group(2, constraints%nh(6), 6)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(6), 6)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(6), 6)= hm_H_index(3, i1)
        constraints%H_Group(5, constraints%nh(6), 6)= hm_H_index(4, i1)
        constraints%H_Group(6, constraints%nh(6), 6)= hm_H_index(5, i1)
        constraints%H_Group(7, constraints%nh(6), 6)= hm_H_index(6, i1)


      else if (hm_duplicate(i1) == 7) then

        constraints%nh(7) = constraints%nh(7) + 1
        constraints%H_Group(1, constraints%nh(7), 7)= i1       
        constraints%H_Group(2, constraints%nh(7), 7)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(7), 7)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(7), 7)= hm_H_index(3, i1)
        constraints%H_Group(5, constraints%nh(7), 7)= hm_H_index(4, i1)
        constraints%H_Group(6, constraints%nh(7), 7)= hm_H_index(5, i1)
        constraints%H_Group(7, constraints%nh(7), 7)= hm_H_index(6, i1)
        constraints%H_Group(8, constraints%nh(7), 7)= hm_H_index(7, i1)


      else if (hm_duplicate(i1) == 8) then

        constraints%nh(8) = constraints%nh(8) + 1
        constraints%H_Group(1, constraints%nh(8), 8)= i1       
        constraints%H_Group(2, constraints%nh(8), 8)= hm_H_index(1, i1)
        constraints%H_Group(3, constraints%nh(8), 8)= hm_H_index(2, i1)
        constraints%H_Group(4, constraints%nh(8), 8)= hm_H_index(3, i1)
        constraints%H_Group(5, constraints%nh(8), 8)= hm_H_index(4, i1)
        constraints%H_Group(6, constraints%nh(8), 8)= hm_H_index(5, i1)
        constraints%H_Group(7, constraints%nh(8), 8)= hm_H_index(6, i1)
        constraints%H_Group(8, constraints%nh(8), 8)= hm_H_index(7, i1)
        constraints%H_Group(9, constraints%nh(8), 8)= hm_H_index(8, i1)

      end if
  
    end do

    return

  end subroutine setup_hbond_group_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_restraints_peer_0
  !> @brief        peer routine of sp_restraints_mod::setup_restraints()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_restraints_peer_0(res_info, sel_info, restraints)

    ! formal arguments
    type(s_res_info),        intent(in)    :: res_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_restraints),      intent(inout) :: restraints

    ! local variables
    type(s_selatoms)         :: selatoms
    integer                  :: i, j, nfunc, ngroup


    ! select atoms
    do i = 1, size(sel_info%groups)
      call select_atom_peer(sel_info%groups(i), selatoms)
      call hm_alloc_selatoms(size(selatoms%idx))
      do j = 1, size(selatoms%idx)
        call hm_set_selatoms(j, i, selatoms%idx(j))
      end do
    end do
    call dealloc_selatoms(selatoms)

    ! copy s_res_info to s_restraints

    restraints%nfunctions = res_info%nfunctions
    restraints%num_groups = size(sel_info%groups)

    nfunc  = restraints%nfunctions
    ngroup = restraints%num_groups

    call alloc_restraints(restraints, RestraintsFunc,  nfunc)
    call alloc_restraints(restraints, RestraintsGroup, ngroup)

    restraints%function     (1:nfunc) = res_info%function     (1:nfunc)
    restraints%exponent     (1:nfunc) = res_info%exponent     (1:nfunc)
    restraints%constant     (1:nfunc) = res_info%constant     (1:nfunc)
    restraints%select_index (1:nfunc) = res_info%select_index (1:nfunc)
    restraints%reference    (1:nfunc) = res_info%reference    (1:nfunc)
    restraints%exponent_dist(1:nfunc) = res_info%exponent_dist(1:nfunc)
    restraints%weight_dist  (1:nfunc) = res_info%weight_dist  (1:nfunc)
    restraints%group       (1:ngroup) = sel_info%groups      (1:ngroup)

    restraints%restraint_flag         = (restraints%nfunctions > 0)

    restraints%max_atoms = 0
    do i = 1, size(sel_info%groups)
      restraints%max_atoms = max(restraints%max_atoms, hm_num_selatoms(i))
    end do

    return

  end subroutine setup_restraints_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_peer_0
  !> @brief        peer routine of sp_enefunc_mod::define_enefunc()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_peer_0(ene_info,   &
                                   par,        &
                                   prmtop,     &
                                   grotop,     &
                                   restraints, &
                                   enefunc)

    ! formal arguments
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_par),          intent(in)    :: par
    type(s_prmtop),       intent(in)    :: prmtop
    type(s_grotop),       intent(in)    :: grotop
    type(s_restraints),   intent(in)    :: restraints
    type(s_enefunc),      intent(inout) :: enefunc


    ! charmm
    !
    if (par%num_bonds > 0) then

      call define_enefunc_charmm_peer_0 (ene_info, par, restraints, enefunc)

    ! amber
    !
    else if (prmtop%num_atoms > 0) then

      call define_enefunc_amber_peer_0  (ene_info, prmtop, restraints, enefunc)

    ! gromacs
    !
    else if (grotop%num_atomtypes > 0) then

      call define_enefunc_gromacs_peer_0(ene_info, grotop, restraints, enefunc)

    end if

    return

  end subroutine define_enefunc_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_charmm_peer_0
  !> @brief        peer routine of 
  !!                            sp_enefunc_charmm_mod::define_enefunc_charmm()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_charmm_peer_0(ene_info,   &
                                          par,        &
                                          restraints, &
                                          enefunc)

    ! formal arguments
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_par),          intent(in)    :: par
    type(s_restraints),   intent(in)    :: restraints
    type(s_enefunc),      intent(inout) :: enefunc


    ! nonbonded
    !
!   call setup_enefunc_nonb_C_peer_0(ene_info, par, enefunc) ! skip


    ! lookup table
    !
    call setup_enefunc_table_peer_0(ene_info, enefunc)


    ! restraint
    !
    call setup_enefunc_restraints_peer_0(restraints, enefunc)

    return

  end subroutine define_enefunc_charmm_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_C_peer_0
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_nonb()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_C_peer_0(ene_info, par, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin
    integer                  :: i, j, nonb_p


    enefunc%num_atom_cls = par%num_atom_cls

    ELECOEF          = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    nonb_p = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nonb_p)

    do i = 1, nonb_p

      lj_coef(1,i) = par%nonb_eps (i)
      lj_coef(2,i) = par%nonb_rmin(i)

      do j = 1, nonb_p

        ! combination rule
        !
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        !
        enefunc%nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
        enefunc%nb14_lj6(i,j)  = 2.0_wp * eps14 * (rmin14 ** 6)
        enefunc%nonb_lj12(i,j) = eps * (rmin ** 12)
        enefunc%nonb_lj6(i,j)  = 2.0_wp * eps * (rmin ** 6)

      end do
    end do

    return

  end subroutine setup_enefunc_nonb_C_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber_peer_0
  !> @brief        peer routine of
  !!                            sp_enefunc_amber_mod::define_enefunc_amber()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber_peer_0(ene_info,   &
                                         prmtop,     &
                                         restraints, &
                                         enefunc)

    ! formal arguments
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_prmtop),       intent(in)    :: prmtop
    type(s_restraints),   intent(in)    :: restraints
    type(s_enefunc),      intent(inout) :: enefunc


    ! nonbonded
    !
!   call setup_enefunc_nonb_A_peer_0(ene_info, prmtop, enefunc) ! skip


    ! lookup table
    !
    call setup_enefunc_table_peer_0(ene_info, enefunc)


    ! restraints
    !
    call setup_enefunc_restraints_peer_0(restraints, enefunc)


    return

  end subroutine define_enefunc_amber_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_A_peer_0
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_nonb()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_A_peer_0(ene_info, prmtop, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, nnonb, ndihetype, ndihe
    integer                  :: nb_par_idx


    ndihetype = size(prmtop%scee_scale_fact)
    call alloc_enefunc(enefunc, EneFuncAMBERScale, ndihetype, 0)

    ELECOEF          = ELECOEF_AMBER

    ndihe = 0
    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_wo_hy(5,i)) =  &
          1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_wo_hy(5,i)) = &
          1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    nnonb = prmtop%num_types

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          enefunc%nb14_lj12(i,j) = 0.0_wp
          enefunc%nb14_lj6 (i,j) = 0.0_wp
          enefunc%nonb_lj12(i,j) = 0.0_wp
          enefunc%nonb_lj6 (i,j) = 0.0_wp
        else
          enefunc%nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          enefunc%nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

      end do
    end do

    enefunc%num_atom_cls = nnonb

    return

  end subroutine setup_enefunc_nonb_A_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs_peer_0
  !> @brief        peer routine of 
  !!                            sp_enefunc_gromacs_mod::define_enefunc_gromacs()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs_peer_0(ene_info,   &
                                           grotop,     &
                                           restraints, &
                                           enefunc)

    ! formal arguments
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_grotop),       intent(in)    :: grotop
    type(s_restraints),   intent(in)    :: restraints
    type(s_enefunc),      intent(inout) :: enefunc


    ! nonbonded
    !
!   call setup_enefunc_nonb_G_peer_0(ene_info, grotop, enefunc) ! skip


    ! lookup table
    !
    call setup_enefunc_table_peer_0(ene_info, enefunc)


    ! restraint
    !
    call setup_enefunc_restraints_peer_0(restraints, enefunc)


    ! restraints (gromacs)
    !
    call setup_enefunc_gro_restraints_peer_0(grotop, enefunc)


    return

  end subroutine define_enefunc_gromacs_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_G_peer_0
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_nonb()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_G_peer_0(ene_info, grotop, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig, ei, ej, si, sj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, ncel, i, j, k, excl_level


    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ELECOEF              = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !

    nnonb = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
      do j = 1, nnonb

        ! combination rule
        !
        if (grotop%num_nbonparms == 0) then

          vi = grotop%atomtypes(i)%v
          vj = grotop%atomtypes(j)%v
          wi = grotop%atomtypes(i)%w
          wj = grotop%atomtypes(j)%w

          if (grotop%defaults%combi_rule == 2) then

            si = vi * 10.0_wp
            sj = vj * 10.0_wp

            ei = wi * JOU2CAL
            ej = wj * JOU2CAL

            sig = (si + sj) * 0.5_wp
            eps = sqrt(ei * ej)

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else

            c6i  = vi * 1000000.0_wp * JOU2CAL
            c6j  = vj * 1000000.0_wp * JOU2CAL

            c12i = wi * 1000000.0_wp * 1000000.0_wp * JOU2CAL
            c12j = wj * 1000000.0_wp * 1000000.0_wp * JOU2CAL

            c6  = sqrt(c6i  * c6j)
            c12 = sqrt(c12i * c12j)

          end if

        else

          vij = 0.0_wp
          wij = 0.0_wp

          do k = 1, grotop%num_nbonparms
            if (grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type2 .or.  &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type2) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (vij == 0.0_wp .or. wij == 0.0_wp) &
       call error_msg('Setup_Enefunc_Nonb_G_Peer_0> combination is not found.')

          if (grotop%defaults%combi_rule == 2) then

            sig = vij * 10.0_wp
            eps = wij * JOU2CAL

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule = 1 or 3

            c6  = vij * 1000000.0_wp * JOU2CAL
            c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL

          end if

        end if

        ! set parameters
        !
        enefunc%nb14_lj12(i,j) = c12
        enefunc%nb14_lj6 (i,j) = c6

        enefunc%nonb_lj12(i,j) = c12
        enefunc%nonb_lj6 (i,j) = c6

      end do
    end do


    ! check # of exclusion level
    !
    if (enefunc%forcefield == ForcefieldGROMARTINI) then

      !TODO

      enefunc%excl_level = -1

      do i = 1, grotop%num_molss
        excl_level = grotop%molss(i)%moltype%exclude_nbon
        if (enefunc%excl_level == -1) then
          enefunc%excl_level = excl_level
        else if (enefunc%excl_level /= excl_level) then
          call error_msg( &
      'Setup_Enefunc_Nonb_G_Peer_0> multiple "exclude_nbon" is not supported.')
        end if
      end do
    end if


    return

  end subroutine setup_enefunc_nonb_G_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_table_peer_0
  !> @brief        peer routine of sp_enefunc_table_mod::setup_enefunc_table()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_table_peer_0(ene_info, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_enefunc),         intent(inout) :: enefunc


    enefunc%table%table_order = ene_info%table_order
    enefunc%table%water_model = ene_info%water_model

    return

  end subroutine setup_enefunc_table_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints_peer_0
  !> @brief        peer routine of 
  !!                       sp_enefunc_restraints_mod::setup_enefunc_restraints()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_restraints_peer_0(restraints, enefunc)

    ! formal arguments
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j
    integer                  :: nconst, nref


    enefunc%restraint = restraints%restraint_flag

    if (.not. restraints%restraint_flag) return

    ! setup group
    !
    call setup_enefunc_rest_group_peer_0(restraints, enefunc)


    ! setup func
    !
    call setup_enefunc_rest_func_peer_0(restraints, enefunc)


    ! setup about domain
    !
    ! ..skip


    if (main_rank) then

      write(MsgOut,'(A)') &
           'Setup_Enefunc_Restraints_Peer_0> Setup restraint groups'

      do i = 1, enefunc%num_restraintgroups

        write(MsgOut,'(A,I5,3A)') &
             ' group = ', i, ', "', trim(restraints%group(i)),'"'
        write(MsgOut,'(A,I5)')    &
             ' # of atoms = ', hm_num_selatoms(i)
        write(MsgOut,'(A)') &
             ' atomlist: '

        do j = 1, hm_num_selatoms(i)

          write(MsgOut,'(i5,$)') hm_selatoms(j,i)
          if (mod(j, 10) == 0 .and. j /= hm_num_selatoms(i)) then
            write(MsgOut,'(A)') ''
          end if

        end do
        write(MsgOut,'(A)') ''

      end do

      write(MsgOut,'(A)') &
           'Setup_Enefunc_Restraints_Peer_0> Setup restraint functions'

      do i = 1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
             ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) then
          write(MsgOut,'(A,4F8.3)') &
               ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        else
          write(MsgOut,'(A,2F8.3)') &
               ' const(Lower, Upper) = ',   enefunc%restraint_const(1:2,i)
          write(MsgOut,'(A,2F8.3)') &
               ' ref  (Lower, Upper) = ',   enefunc%restraint_ref(1:2,i)
        end if

        write(MsgOut,'(A,I5)') ' exponend of function = ',                &
             enefunc%restraint_exponent_func(i)

        ! summary for distance restraint
        !
        if (enefunc%restraint_kind(i) == RestraintsFuncDIST .or.          &
            enefunc%restraint_kind(i) == RestraintsFuncDISTCOM) then

          call error_msg('Setup_Enefunc_Restraints_Peer_0> '//&
               'DIST restraint is not supported in parallel I/O')

        ! summary for angle restraint
        !
        else if ( &
            enefunc%restraint_kind(i) == RestraintsFuncANGLE .or.    &
            enefunc%restraint_kind(i) == RestraintsFuncANGLECOM) then

          call error_msg('Setup_Enefunc_Restraints_Peer_0> '//&
               'ANGLE restraint is not supported in parallel I/O')

        ! summary for dihedral angle restraint
        !
        else if ( &
            enefunc%restraint_kind(i) == RestraintsFuncDIHED  .or.   &
            enefunc%restraint_kind(i) == RestraintsFuncDIHEDCOM) then

          call error_msg('Setup_Enefunc_Restraints_Peer_0> '//&
               'DIHED restraint is not supported in parallel I/O')

        ! summary for pc restraint
        !
        else if ( &
            enefunc%restraint_kind(i) == RestraintsFuncPC  .or.   &
            enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then

          call error_msg('Setup_Enefunc_Restraints_Peer_0> '//&
               'PC restraint is not supported in parallel I/O')

        ! summary for positinal and RMSD restraint
        !
        else 

          write(MsgOut,'(A,I5)') ' # of groups  = ',enefunc%restraint_funcgrp(i)
          write(MsgOut,'(" grouplist: ",$)') 
          do j = 1, enefunc%restraint_funcgrp(i) 
            write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
            if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
                 write(MsgOut,'(A)') ''
          end do
          write(MsgOut,'(A)') ''

        end if
        write(MsgOut,'(A)') ''

      end do

    end if

    return

  end subroutine setup_enefunc_restraints_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_group_peer_0
  !> @brief        peer routine of 
  !!                       sp_enefunc_restraints_mod::setup_enefunc_rest_group()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_group_peer_0(restraints, enefunc)

    ! formal arguments
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc


    enefunc%num_restraintgroups = restraints%num_groups

    return

  end subroutine setup_enefunc_rest_group_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_enefunc_rest_func_peer_0
  !> @brief        peer routine of 
  !!                        sp_enefunc_restraints_mod::setup_enefunc_rest_func()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_func_peer_0(restraints, enefunc) 

    ! formal arguments
    type(s_restraints),      intent(in)     :: restraints
    type(s_enefunc),         intent(inout)  :: enefunc

    ! local variables
    integer                  :: i, ndata
    integer                  :: num_funcs, max_grp
    integer,     allocatable :: num_indexes(:)


    ! get total number of groups in index
    !
    enefunc%num_restraintfuncs = restraints%nfunctions
    num_funcs = enefunc%num_restraintfuncs
    allocate(num_indexes(1:num_funcs))

    do i = 1, num_funcs
      ndata = split_num(restraints%select_index(i))

      select case(restraints%function(i))
      case(RestraintsFuncDIST:RestraintsFuncDISTCOM)
        if (ndata <  2 .or. mod(ndata,2) /= 0) &
          call error_msg( &
            'Setup_Enefunc_Rest_Func_Peer_0> Error in index for DIST')

      case(RestraintsFuncANGLE:RestraintsFuncANGLECOM)
        if (ndata /= 3) &
          call error_msg( &
            'Setup_Enefunc_Rest_Func_Peer_0> Error in index for ANGLE')

      case(RestraintsFuncDIHED:RestraintsFuncDIHEDCOM)
        if (ndata /= 4) &
          call error_msg( &
            'Setup_Enefunc_Rest_Func_Peer_0> Error in index for DIHE')

      case default
        if (ndata /= 1) & 
          call error_msg( &
            'Setup_Enefunc_Rest_Func_Peer_0> Error in index for POSI or RMSD')

      end select

      num_indexes(i) = ndata
    end do

    ! get maximum number of groups among all 'index's
    !
    max_grp = maxval(num_indexes(1:num_funcs))
    enefunc%max_restraint_numgrps = max_grp

    call alloc_enefunc(enefunc, EneFuncReff, num_funcs, max_grp)


    ! setup parameters
    !
    do i = 1, num_funcs

      enefunc%restraint_funcgrp(i) = num_indexes(i)
      enefunc%restraint_kind(i)    = restraints%function(i)

      ! setup group list
      !
      call setup_restraints_grouplist(i, restraints, enefunc)

      ! setup exponent_func
      !
      call setup_restraints_exponent_func(i, restraints, enefunc)

      ! setup force constants
      !
      call setup_restraints_constants(i, restraints, enefunc)

      ! setup reference   (except for POSI)
      !
      call setup_restraints_reference(i, restraints, enefunc)

      ! setup exponent_dist (for DIST only)
      !
      call setup_restraints_exponent_dist(i, restraints, enefunc)

      ! setup weight_dist   (for DIST only)
      !
      call setup_restraints_weight_dist(i, restraints, enefunc)


    end do

    ! setup reference coordinates (for POSI and RMSD)
    !
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.    &
          enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncPC .or.      &
          enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then

        call setup_enefunc_rest_refcoord_peer_0(enefunc)

        if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) &
          enefunc%restraint_posi = .true.
        if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or.  &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM)   &
          enefunc%restraint_rmsd = .true.

      end if

    end do

    ! deallocate local array
    !
    deallocate(num_indexes)

    return

  end subroutine setup_enefunc_rest_func_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_enefunc_gro_restraints_peer_0
  !> @brief        peer routine of 
  !!                      sp_enefunc_gromacs_mod::setup_enefunc_gro_restraints()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gro_restraints_peer_0(grotop, enefunc) 

    ! formal arguments
    type(s_grotop),          intent(in)     :: grotop
    type(s_enefunc),         intent(inout)  :: enefunc

    ! local variables
    type(s_enefunc)          :: ef0
    type(s_selatoms)         :: selatoms
    real(wp)                 :: kx, ky, kz
    integer                  :: nposres, npr_atom, max_pr_atom, ioffset
    integer                  :: i, j, k, n, n2, n3, group0, func0, istart, iend
    logical                  :: lexist

    type(s_grotop_mol), pointer :: gromol


    ! count # of position restraints
    !
    nposres     = 0
    max_pr_atom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = grotop%molss(i)%count * gromol%num_posress
      max_pr_atom = max(max_pr_atom,npr_atom)
    end do

    if (nposres == 0) &
      return

    enefunc%restraint = .true.


    ! setup EneFuncRefc
    !
    do i = 1, hm_num_atoms
      call hm_set_atom_refcoord(1, i, hm_atom_coord(1, i))
      call hm_set_atom_refcoord(2, i, hm_atom_coord(2, i))
      call hm_set_atom_refcoord(3, i, hm_atom_coord(3, i))
    end do


    ! setup EneFuncRefg
    !
    group0 = size(enefunc%restraint_numatoms)

    call alloc_selatoms(selatoms, max_pr_atom)

    ioffset = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      npr_atom = 0
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_posress
          npr_atom = npr_atom + 1
          selatoms%idx(npr_atom) = gromol%posress(k)%atom_idx + ioffset

          if (k > 1 .and. main_rank) then
            if (gromol%posress(k-1)%kx /= gromol%posress(k)%kx .or. &
                gromol%posress(k-1)%ky /= gromol%posress(k)%ky .or. &
                gromol%posress(k-1)%kz /= gromol%posress(k)%kz) then
              write(MsgOut,'(a)')'Setup_Enefunc_Gro_Restraints_Peer_0> WARNING:'
              write(MsgOut,'(a,a,a)') &
  '   different restraint constant between foreach atoms is not supported. [', &
  trim(grotop%molss(i)%moltype%name), ']'
              write(MsgOut,'(a)') ' '
            end if
          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

      call hm_alloc_selatoms(npr_atom)
      do j = 1, npr_atom
        call hm_set_selatoms(j, i, selatoms%idx(j))
      end do

    end do

    call dealloc_selatoms(selatoms)


    ! setup EneFuncReff
    !
    if (allocated(enefunc%restraint_kind)) then

      n  = size(enefunc%restraint_kind)
      n2 = size(enefunc%restraint_grouplist(:,1))
      n3 = max(int(n2/2),1)
      call alloc_enefunc(ef0, EneFuncReff, n, n2)
      ef0%restraint_kind          (1:n) = enefunc%restraint_kind          (1:n)
      ef0%restraint_grouplist(1:n2,1:n) = enefunc%restraint_grouplist(1:n2,1:n)
      ef0%restraint_const    (1:4, 1:n) = enefunc%restraint_const    (1:4, 1:n)
      ef0%restraint_ref      (1:2, 1:n) = enefunc%restraint_ref      (1:2, 1:n)
      ef0%restraint_funcgrp       (1:n) = enefunc%restraint_funcgrp       (1:n)
      ef0%restraint_exponent_func (1:n) = enefunc%restraint_exponent_func (1:n)
      ef0%restraint_exponent_dist (1:n3,1:n) &
                                   = enefunc%restraint_exponent_dist (1:n3,1:n)
      ef0%restraint_weight_dist   (1:n3,1:n) &
                                   = enefunc%restraint_weight_dist   (1:n3,1:n)

      call alloc_enefunc(enefunc, EneFuncReff, n+nposres, max(n2,1))
      enefunc%restraint_kind          (1:n) = ef0%restraint_kind          (1:n)
      enefunc%restraint_grouplist(1:n2,1:n) = ef0%restraint_grouplist(1:n2,1:n)
      enefunc%restraint_const    (1:4, 1:n) = ef0%restraint_const    (1:4, 1:n)
      enefunc%restraint_ref      (1:2, 1:n) = ef0%restraint_ref      (1:2, 1:n)
      enefunc%restraint_funcgrp       (1:n) = ef0%restraint_funcgrp       (1:n)
      enefunc%restraint_exponent_func (1:n) = ef0%restraint_exponent_func (1:n)
      enefunc%restraint_exponent_dist (1:n3,1:n) &
                                       = ef0%restraint_exponent_dist (1:n3,1:n)
      enefunc%restraint_weight_dist   (1:n3,1:n) &
                                       = ef0%restraint_weight_dist   (1:n3,1:n)

      call dealloc_enefunc(ef0, EneFuncReff)

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncReff, nposres, 1)

    end if


    func0 = n

    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress == 0) &
        cycle

      nposres = nposres + 1

      kx = gromol%posress(1)%kx * 0.01_wp * JOU2CAL * 0.5_wp
      ky = gromol%posress(1)%ky * 0.01_wp * JOU2CAL * 0.5_wp
      kz = gromol%posress(1)%kz * 0.01_wp * JOU2CAL * 0.5_wp

      enefunc%restraint_kind         (func0+nposres) = RestraintsFuncPOSI
      enefunc%restraint_funcgrp      (func0+nposres) = 1
      enefunc%restraint_grouplist  (1,func0+nposres) = group0 + nposres
      enefunc%restraint_const      (1,func0+nposres) = 1.0_wp
      enefunc%restraint_const      (2,func0+nposres) = kx
      enefunc%restraint_const      (3,func0+nposres) = ky
      enefunc%restraint_const      (4,func0+nposres) = kz
      enefunc%restraint_exponent_func(func0+nposres) = 2

    end do

    enefunc%num_restraintgroups = group0+nposres
    enefunc%num_restraintfuncs  = func0+nposres


    ! summary of setup enefunc_gro_restraints
    !
    if (main_rank) then

      write(MsgOut,'(A)')&
           'Setup_Enefunc_Gro_Restraints_Peer_0> Setup restraint functions'

      do i = func0+1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        write(MsgOut,'(A,4F8.3)') &
             ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        write(MsgOut,'(A,I5)') ' exponend of function = ', &
                                          enefunc%restraint_exponent_func(i)

        write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%restraint_funcgrp(i)
        write(MsgOut,'(" grouplist: ",$)') 
        do j = 1, enefunc%restraint_funcgrp(i) 
          write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
          if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
            write(MsgOut,'(A)') ''
        end do
        write(MsgOut,'(A)') ''

      end do

      write(MsgOut,*) ''

    end if

    return

  end subroutine setup_enefunc_gro_restraints_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_enefunc_rest_refcoord_peer_0
  !> @brief        peer routine of 
  !!                    sp_enefunc_restraints_mod::setup_enefunc_rest_refcoord()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_refcoord_peer_0(enefunc) 

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc


    if (.not. hm_atom_ref_set_flag) &
      call error_msg('Setup_Enefunc_Refcoord_Peer_0> Refcoord is not set')

    return

  end subroutine setup_enefunc_rest_refcoord_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics_peer_0
  !> @brief        peer routine of sp_dynamics_mod::setup_dynamics()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics_peer_0(dyn_info, bound_info, res_info, dynamics)

    ! formal arguments
    type(s_dyn_info),        intent(in)    :: dyn_info
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_res_info),        intent(in)    :: res_info
    type(s_dynamics),        intent(inout) :: dynamics

    integer                                :: iseed


    ! initialize variables in dynamics
    !
    call init_dynamics(dynamics)


    ! setup variables
    !
    iseed = dyn_info%iseed
    if (dyn_info%iseed == -1) then
      call error_msg('Setup_Dynamics_Peer_0> prst_setup is not '//&
                     'allowed iseed=-1')
    endif
    dynamics%iseed = iseed

    ! setup random system
    !
    call random_init(dynamics%iseed)
    call random_push_stock

    return

  end subroutine setup_dynamics_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints_peer_0
  !> @brief        peer routine of sp_constraints_mod::setup_constraints()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints_peer_0(cons_info, par, prmtop, grotop, &
                                      enefunc, constraints)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    integer                  :: iw1, iw2


    constraints%rigid_bond      = cons_info%rigid_bond
    constraints%fast_bond       = cons_info%fast_bond
    constraints%fast_water      = cons_info%fast_water
    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_order     = cons_info%lincs_order
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%water_model     = cons_info%water_model

    if (constraints%rigid_bond) then

      ! setup SETTLE
      !
      if (constraints%fast_water) then

        if (constraints%tip4) then
          call setup_fast_water_tip4_peer_0(par, prmtop, grotop, &
                                            enefunc, constraints)
        else
          call setup_fast_water_peer_0(par, prmtop, grotop,      &
                                       enefunc, constraints)
        end if

      end if


      ! setup SHAKE and RATTLE
      !
      call setup_rigid_bond_peer_0(par, enefunc, constraints)

    end if

    return

  end subroutine setup_constraints_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_peer_tip4_0
  !> @brief        peer routine of sp_constraints_mod::setup_fast_water_tip4()
  !! @authors      NT, JJ
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_tip4_peer_0(par, prmtop, grotop,  enefunc, &
                                          constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, i1, i2, list(4), ioffset
    character(6)             :: ci1, ci2, ci3, ci4

    type(s_grotop_mol), pointer :: gromol


    ! mass
    !
    constraints%water_massO = enefunc%table%mass_O
    constraints%water_massH = enefunc%table%mass_H

    ! min distance
    !

    list(1) = hm_water_list(1,1)
    list(2) = hm_water_list(2,1)
    list(3) = hm_water_list(3,1)
    list(4) = hm_water_list(4,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))
      ci4 = hm_atom_cls_name(list(4))

      do i = 1, par%num_bonds
        if ((ci1 == par%bond_atom_cls(1,i) .and. &
             ci2 == par%bond_atom_cls(2,i)) .or. &
            (ci1 == par%bond_atom_cls(2,i) .and. &
             ci2 == par%bond_atom_cls(1,i))) then
          constraints%water_rOH = par%bond_dist_min(i)
          exit
        end if
      end do

      do i = 1, par%num_bonds
        if (ci2 == par%bond_atom_cls(1,i) .and. &
            ci3 == par%bond_atom_cls(2,i)) then
          constraints%water_rHH = par%bond_dist_min(i)
          exit
        end if
      end do

      do i = 1, par%num_bonds
        if ((ci1 == par%bond_atom_cls(1,i) .and. &
             ci4 == par%bond_atom_cls(2,i)) .or. &
            (ci4 == par%bond_atom_cls(2,i) .and. &
             ci1 == par%bond_atom_cls(1,i))) then
          constraints%water_rOD = par%bond_dist_min(i)
          exit
        end if
      end do

    ! amber
    else if (prmtop%num_atoms > 0) then

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(2) == i2 .or.  &
            list(2) == i1 .and. list(1) == i2) then

          constraints%water_rOH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(2) == i1 .and. list(3) == i2 .or. &
            list(3) == i1 .and. list(2) == i2) then

          constraints%water_rHH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_mbonda

        i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
        i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(4) == i2 .or.  &
            list(4) == i1 .and. list(1) == i2) then

          constraints%water_rOD = &
               prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
          exit

        end if

      end do

    end if

    return

  end subroutine setup_fast_water_tip4_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_peer_0
  !> @brief        peer routine of sp_constraints_mod::setup_fast_water()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_peer_0(par, prmtop, grotop,  enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, j, k, i1, i2, list(3), ioffset
    character(6)             :: ci1, ci2, ci3

    type(s_grotop_mol), pointer :: gromol


    ! mass
    !
    constraints%water_massO = enefunc%table%mass_O
    constraints%water_massH = enefunc%table%mass_H

    ! min distance
    !

    list(1) = hm_water_list(1,1)
    list(2) = hm_water_list(2,1)
    list(3) = hm_water_list(3,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))

      do i = 1, par%num_bonds
        if ((ci1 == par%bond_atom_cls(1,i) .and. &
             ci2 == par%bond_atom_cls(2,i)) .or. &
             (ci1 == par%bond_atom_cls(2,i) .and. &
             ci2 == par%bond_atom_cls(1,i))) then
          constraints%water_rOH = par%bond_dist_min(i)
          exit
        end if
      end do

      do i = 1, par%num_bonds
        if (ci2 == par%bond_atom_cls(1,i) .and. &
             ci3 == par%bond_atom_cls(2,i)) then
          constraints%water_rHH = par%bond_dist_min(i)
          exit
        end if
      end do

    ! amber
    else if (prmtop%num_atoms > 0) then

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(1) == i1 .and. list(2) == i2 .or.  &
            list(2) == i1 .and. list(1) == i2) then

          constraints%water_rOH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

      do i = 1, prmtop%num_bondh

        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

        if (list(2) == i1 .and. list(3) == i2 .or. &
            list(3) == i1 .and. list(2) == i2) then

          constraints%water_rHH = &
               prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
          exit

        end if

      end do

    ! gromacs
    else if (grotop%num_atomtypes > 0) then

      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(1) == i1 .and. list(2) == i2 .or.  &
                  list(2) == i1 .and. list(1) == i2) then

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_wp
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_wp
          goto 1

        end if

      end do

1     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        if (gromol%settles%func == 0) then

          do j = 1, grotop%molss(i)%count

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1 + ioffset
              i2 = gromol%bonds(k)%atom_idx2 + ioffset

              if (list(2) == i1 .and. list(3) == i2 .or.  &
                  list(3) == i1 .and. list(2) == i2) then

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_wp
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rHH = gromol%settles%dhh * 10.0_wp
          goto 2
          
        end if

      end do
2     continue

    end if

    return

  end subroutine setup_fast_water_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigid_bond_peer_0
  !> @brief        peer routine of sp_constraints_mod::setup_rigid_bond()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigid_bond_peer_0(par, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints


    ! nothing to do

    return

  end subroutine setup_rigid_bond_peer_0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_md_peer
  !> @brief        peer routine of sp_setup_spdyn_mod::setup_md()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_md_peer(local_my_rank, &
                           ctrl_data,     &
                           par,           &
                           prmtop,        &
                           grotop,        &
                           domain_index,  &
                           boundary,      &
                           enefunc,       &
                           constraints,   &
                           domain)

    ! formal arguments
    integer,                 intent(in)    :: local_my_rank
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain_index),    intent(in)    :: domain_index
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! set parameters for domain 
    !
    call setup_domain_peer(ctrl_data%ene_info,  &
                           ctrl_data%cons_info, &
                           local_my_rank, domain_index, boundary, &
                           enefunc, constraints, domain)

    ! setup enefunc in each domain
    !
    call define_enefunc_peer(par, prmtop, grotop, ctrl_data%ene_info, &
                             domain_index, constraints, domain, enefunc)

    return

  end subroutine setup_md_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_peer
  !> @brief        peer routine of sp_domain_mod::setup_domain()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_peer(ene_info,      &
                               con_info,      &
                               local_my_rank, &
                               domain_index,  &
                               boundary,      &
                               enefunc,       &
                               constraints,   &
                               domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    integer,                 intent(in)    :: local_my_rank
    type(s_domain_index),    intent(in)    :: domain_index
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all
    real(wp)                 :: mass_upper_bound


    ! initialize structure informations
    !
    call init_domain(domain)

    domain%num_atom_all = hm_num_atoms

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model

    mass_upper_bound = con_info%hydrogen_mass_upper_bound

    ! assign the rank of each dimension from my_rank
    !

    my_world_rank = local_my_rank
    my_city_rank  = local_my_rank

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

    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all,  1, 1)
    call alloc_id_g2l(domain, hm_num_atoms)

    if (constraints%rigid_bond) then

      if (constraints%tip4) then
        call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
      else
        call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
      end if
      call setup_atom_by_HBond_peer(mass_upper_bound, domain_index, &
                                    boundary, constraints, domain)

    else
      if (enefunc%table%tip4) then
        call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
      else
        call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
      end if

      call setup_atom_by_table_peer(domain_index, boundary, enefunc, domain)

    end if


    my_world_rank = local_my_rank
    my_city_rank  = local_my_rank

    call setup_domain_interaction(boundary, domain)

    return

  end subroutine setup_domain_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_HBond_peer
  !> @brief        peer routine of sp_domain_mod::setup_atom_by_HBond_peer()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_HBond_peer(mass_upper_bound, &
                                      domain_index,     &
                                      boundary,         &
                                      constraints,      &
                                      domain)

    ! formal arguments
    real(wp),                     intent(in)    :: mass_upper_bound
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_boundary),     target, intent(in)    :: boundary
    type(s_constraints),  target, intent(inout) :: constraints
    type(s_domain),       target, intent(inout) :: domain

    ! local variable
    real(wp)                      :: x_shift, y_shift, z_shift
    real(wp)                      :: move(3), origin(3), mass
    integer                       :: i, j, icx, icy, icz, icel
    integer                       :: isolute, iwater, ih1, ih2, id
    integer                       :: iisolute, iiwater
    integer                       :: icel_local
    integer                       :: ncel, ncel_local
    character(6)                  :: ci1, ci2
    logical                       :: cl1, cl2, mi1, mi2
    logical                       :: tip4

    real(wp),             pointer :: bsize_x, bsize_y, bsize_z
    real(wp),             pointer :: csize_x, csize_y, csize_z
    integer,              pointer :: iwater_list(:), isolute_list(:)
    integer,              pointer :: ncel_x, ncel_y, ncel_z
    integer,              pointer :: cell_g2l(:), cell_g2b(:)
    integer,              pointer :: natom(:), nsolute(:), nwater(:)
    integer,              pointer :: solute_list(:,:), water_list(:,:,:)
    integer,              pointer :: No_HGr(:), HGr_local(:,:)
    integer,              pointer :: HGr_bond_list(:,:,:,:)
    integer,              pointer :: patm, psol, pwat, pnoh, phgl


    ncel          = domain%num_cell_local + domain%num_cell_boundary
    ncel_local    = domain%num_cell_local

    call alloc_constraints(constraints, ConstraintsDomainBond, ncel, &
                           constraints%connect)

    iwater_list   => domain_index%iwater_list
    isolute_list  => domain_index%isolute_list

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

    origin(1)     = boundary%origin_x
    origin(2)     = boundary%origin_y
    origin(3)     = boundary%origin_z

    tip4          = constraints%tip4

    natom      (                       1:ncel) = 0
    nsolute    (                       1:ncel) = 0
    nwater     (                       1:ncel) = 0
    solute_list(          1:MaxAtom,   1:ncel) = 0
    water_list (1:3,      1:MaxWater,  1:ncel) = 0
    No_HGr     (                       1:ncel) = 0
    HGr_local  (1:constraints%connect, 1:ncel) = 0

    ! solute atoms (not bonded to hydrogen) in each domain
    !
    do iisolute = 1, domain_index%nsolute_list

      isolute = isolute_list(iisolute)

      i = hm_solute_list(isolute)

      mass = hm_mass(i)
      if (mass > EPS .and. mass < mass_upper_bound) then
        mi1 = .true.
      else 
        mi1 = .false.
      end if
      cl1 = .false.
      ci1 = hm_atom_name(i)
      if (ci1(1:1) == 'H' .or. ci1(1:1) == 'h' .or. &
          ci1(1:1) == 'D' .or. ci1(1:1) == 'd') then
        cl1 = .true.
      end if
      ci1 = hm_atom_cls_name(i)
      if (ci1(1:1) == 'H' .or. ci1(1:1) == 'h' .or. &
          ci1(1:1) == 'D' .or. ci1(1:1) == 'd') then
        cl1 = .true.
      end if
      
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1)
      end if

      if (hm_duplicate(i) == 0 .and. .not. cl1) then

        !coordinate shifted against the origin
        !
        x_shift = hm_atom_coord(1, i) - boundary%origin_x
        y_shift = hm_atom_coord(2, i) - boundary%origin_y
        z_shift = hm_atom_coord(3, i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
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

          patm => natom  (icel_local)
          psol => nsolute(icel_local)
          pnoh => No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain_peer(move, origin, i, &
                                       domain, icel_local, patm)

        ! atoms in the boundary
        !
        else if (cell_g2b(icel) /= 0) then
  
          ! local cell index
          !
          icel_local = cell_g2b(icel) + ncel_local

          patm => natom  (icel_local)
          psol => nsolute(icel_local)
          pnoh => No_HGr (icel_local)

          ! local_count : total number of atoms in each cell
          !
          patm = patm + 1
          psol = psol + 1
          pnoh = pnoh + 1

          solute_list(psol,icel_local) = patm

          call molecule_to_domain_peer(move, origin, i, &
                                       domain, icel_local, patm)

        end if 

      end if

    end do

    ! Hydrogen bonding group
    !
    do isolute = 1, constraints%connect

      do j = 1, constraints%nh(isolute)

        i = constraints%H_Group(1, j, isolute)
        x_shift = hm_atom_coord(1, i) - boundary%origin_x
        y_shift = hm_atom_coord(2, i) - boundary%origin_y
        z_shift = hm_atom_coord(3, i) - boundary%origin_z

        move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
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

          patm => natom  (icel_local)
          psol => nsolute(icel_local)
          phgl => HGr_local(isolute,icel_local)
  
          ! nonhydrogen atom
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm 

          call molecule_to_domain_peer(move, origin, i, &
                                       domain, icel_local, patm)

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1, j, isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm 

            call molecule_to_domain_peer(move, origin, i, &
                                         domain, icel_local, patm)

          end do

        else if (cell_g2b(icel) /= 0) then

          icel_local = cell_g2b(icel) + ncel_local

          patm => natom  (icel_local)
          psol => nsolute(icel_local)
          phgl => HGr_local(isolute,icel_local)

          ! hydrogen atoms
          !
          patm = patm + 1
          psol = psol + 1
          phgl = phgl + 1

          solute_list(psol,icel_local) = patm
          HGr_bond_list(1,phgl,isolute,icel_local) = patm 

          call molecule_to_domain_peer(move, origin, i, &
                                       domain, icel_local, patm)

          do ih1 = 1, isolute

            i = constraints%H_Group(ih1+1,j,isolute)

            patm = patm + 1
            psol = psol + 1

            solute_list(psol,icel_local) = patm
            HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm 

            call molecule_to_domain_peer(move, origin, i, &
                                         domain, icel_local, patm)

          end do

        end if

      end do

    end do

    ! water atoms in each domain
    !
    do iiwater = 1, domain_index%nwater_list

      iwater = iwater_list(iiwater)

      i   = hm_water_list(1,iwater)
      ih1 = hm_water_list(2,iwater)
      ih2 = hm_water_list(3,iwater)
      if (tip4) id = hm_water_list(4,iwater)

      !coordinate shifted against the origin
      !
      x_shift = hm_atom_coord(1, i) - boundary%origin_x
      y_shift = hm_atom_coord(2, i) - boundary%origin_y
      z_shift = hm_atom_coord(3, i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
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

        patm => natom (icel_local)
        pwat => nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih1, &
                                     domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih2, &
                                     domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_peer(move, origin, id, &
                                       domain, icel_local, patm)
        end if

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm => natom (icel_local)
        pwat => nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih1, &
                                     domain, icel_local, patm)


        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3,pwat,icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih2, &
                                     domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1
          water_list(4,pwat,icel_local) = patm

          call molecule_to_domain_peer(move, origin, id, &
                                       domain, icel_local, patm)

        end if

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_HBond_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_table_peer
  !> @brief        peer routine of sp_domain_mod::setup_atom_table()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_table_peer(domain_index, &
                                      boundary,     &
                                      enefunc,      &
                                      domain)

    ! formal arguments
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_boundary),     target, intent(in)    :: boundary
    type(s_enefunc),      target, intent(in)    :: enefunc 
    type(s_domain),       target, intent(inout) :: domain

    ! local variables
    real(wp)                      :: x_shift, y_shift, z_shift
    real(wp)                      :: move(3), origin(3)
    integer                       :: i, ih1, ih2, id
    integer                       :: icx, icy, icz, icel, icel_local
    integer                       :: isolute, iwater
    integer                       :: iisolute, iiwater
    integer                       :: ncel, ncel_local
    logical                       :: tip4

    real(wp),             pointer :: bsize_x, bsize_y, bsize_z
    real(wp),             pointer :: csize_x, csize_y, csize_z
    integer,              pointer :: iwater_list(:), isolute_list(:)
    integer,              pointer :: ncel_x, ncel_y, ncel_z
    integer,              pointer :: cell_g2l(:), cell_g2b(:)
    integer,              pointer :: natom(:), nsolute(:), nwater(:)
    integer,              pointer :: solute_list(:,:), water_list(:,:,:)
    integer,              pointer :: patm, psol, pwat


    iwater_list  => domain_index%iwater_list
    isolute_list => domain_index%isolute_list

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

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local

    tip4         = enefunc%table%tip4

    natom      (                 1:ncel) = 0
    nsolute    (                 1:ncel) = 0
    nwater     (                 1:ncel) = 0
    solute_list(     1:MaxAtom,  1:ncel) = 0
    water_list (1:3, 1:MaxWater, 1:ncel) = 0


    ! solute atoms in each domain
    !
    do iisolute = 1, domain_index%nsolute_list

      isolute = isolute_list(iisolute)

      i = hm_solute_list(isolute)

      ! coordinate shifted against the origin
      !
      x_shift = hm_atom_coord(1,i) - boundary%origin_x
      y_shift = hm_atom_coord(2,i) - boundary%origin_y
      z_shift = hm_atom_coord(3,i) - boundary%origin_z

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      ! assign which cell
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

        patm => natom  (icel_local)
        psol => nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1

        solute_list(psol, icel_local) = patm
        solute_list(patm, icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm => natom  (icel_local)
        psol => nsolute(icel_local)

        ! local_count : total number of atoms in each cell
        !
        patm = patm + 1
        psol = psol + 1

        solute_list(psol, icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

      end if

    end do


    ! water atoms in each domain
    !
    do iiwater = 1, domain_index%nwater_list

      iwater = iwater_list(iiwater)

      i   = hm_water_list(1, iwater)
      ih1 = hm_water_list(2, iwater)
      ih2 = hm_water_list(3, iwater)
      if (tip4) id = hm_water_list(4, iwater)

      ! coordinate shifted against the origin
      !
      x_shift = hm_atom_coord(1, i) - boundary%origin_x
      y_shift = hm_atom_coord(2, i) - boundary%origin_y
      z_shift = hm_atom_coord(3, i) - boundary%origin_z

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
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

        patm => natom (icel_local)
        pwat => nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih1, &
                                     domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih2, &
                                     domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1 
          water_list(4, pwat, icel_local) = patm
          call molecule_to_domain_peer(move, origin, id, &
                                       domain, icel_local, patm)
        end if

      ! atoms in the boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        patm => natom (icel_local)
        pwat => nwater(icel_local)

        pwat = pwat + 1

        ! oxygen atoms
        !
        patm = patm + 1
        water_list(1, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, patm)

        ! first hydrogen atoms
        !
        patm = patm + 1
        water_list(2, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih1, &
                                     domain, icel_local, patm)

        ! second hydrogen atoms
        !
        patm = patm + 1
        water_list(3, pwat, icel_local) = patm

        call molecule_to_domain_peer(move, origin, ih2, &
                                     domain, icel_local, patm)

        ! dummy atoms
        !
        if (tip4) then
          patm = patm + 1 
          water_list(4, pwat, icel_local) = patm
          call molecule_to_domain_peer(move, origin, id, &
                                       domain, icel_local, patm)
        end if

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_by_table_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_peer
  !> @brief        peer routine of sp_domain_mod::setup_atom()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_peer(domain_index, &
                             boundary,     &
                             domain)

    ! formal arguments
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_boundary),     target, intent(in)    :: boundary
    type(s_domain),       target, intent(inout) :: domain

    ! local variables
    real(wp)                      :: x_shift, y_shift, z_shift
    real(wp)                      :: move(3), origin(3)
    integer                       :: i, iatom
    integer                       :: icx, icy, icz, icel, icel_local
    integer                       :: ncel, ncel_local

    real(wp),             pointer :: bsize_x, bsize_y, bsize_z
    real(wp),             pointer :: csize_x, csize_y, csize_z
    integer,              pointer :: iatom_list(:)
    integer,              pointer :: ncel_x, ncel_y, ncel_z
    integer,              pointer :: cell_g2l(:), cell_g2b(:)
    integer,              pointer :: natom(:)


    iatom_list   => domain_index%iatom_list

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

    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    ncel_local   = domain%num_cell_local

    natom      (                 1:ncel) = 0


    do iatom = 1, domain_index%natom_list

      i = iatom_list(iatom)

      !coordinate shifted against the origin
      !
      x_shift = hm_atom_coord(1,i) - boundary%origin_x
      y_shift = hm_atom_coord(2,i) - boundary%origin_y
      z_shift = hm_atom_coord(3,i) - boundary%origin_z

      !coordinate shifted to the first quadrant and set into the boundary box
      !
      move(1) = bsize_x*0.5_wp - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wp - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wp - bsize_z*anint(z_shift/bsize_z)
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

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, natom(icel_local))


      ! atoms in a boundary
      !
      else if (cell_g2b(icel) /= 0) then

        ! local cell index
        !
        icel_local = cell_g2b(icel) + ncel_local

        natom(icel_local) = natom(icel_local) + 1

        call molecule_to_domain_peer(move, origin, i, &
                                     domain, icel_local, natom(icel_local))

      end if

    end do

    domain%num_atom_t0(1:ncel) = natom(1:ncel)

    return

  end subroutine setup_atom_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain_peer
  !> @brief        peer routine of sp_domain_mod::molecule_to_domain_peer()
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain_peer(move, origin, iatom, &
                                     domain, icel, icel_atom)

    ! formal arguments
    real(wp),                 intent(in)    :: move(3)
    real(wp),                 intent(in)    :: origin(3)
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    ! local variables
    real(wp),         pointer :: coord_local (:,:,:)
    real(wp),         pointer :: ref_local   (:,:,:)
    real(wp),         pointer :: vel_local   (:,:,:)
    real(wp),         pointer :: charge_local(:,:)
    real(wp),         pointer :: mass_local  (:,:)
    real(wp),         pointer :: trans       (:,:,:)
    integer,          pointer :: class_local (:,:)
    integer,          pointer :: id_l2g      (:,:)
    integer,          pointer :: id_g2l      (:,:)
    

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

    coord_local(1,icel_atom,icel) = hm_atom_coord(1,iatom)
    coord_local(2,icel_atom,icel) = hm_atom_coord(2,iatom)
    coord_local(3,icel_atom,icel) = hm_atom_coord(3,iatom)
    ref_local  (1,icel_atom,icel) = hm_atom_coord(1,iatom)
    ref_local  (2,icel_atom,icel) = hm_atom_coord(2,iatom)
    ref_local  (3,icel_atom,icel) = hm_atom_coord(3,iatom)
    vel_local  (1,icel_atom,icel) = hm_atom_velocity(1,iatom)
    vel_local  (2,icel_atom,icel) = hm_atom_velocity(2,iatom)
    vel_local  (3,icel_atom,icel) = hm_atom_velocity(3,iatom)
    charge_local (icel_atom,icel) = hm_charge(iatom)
    mass_local   (icel_atom,icel) = hm_mass(iatom)
    class_local  (icel_atom,icel) = hm_atom_cls_no(iatom)
    trans      (1,icel_atom,icel) = move(1) - origin(1)
    trans      (2,icel_atom,icel) = move(2) - origin(2)
    trans      (3,icel_atom,icel) = move(3) - origin(3)

    return

  end subroutine molecule_to_domain_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_peer
  !> @brief        peer routine of sp_enefunc_mod::define_enefunc()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_peer(par,          &
                                 prmtop,       &
                                 grotop,       &
                                 ene_info,     &
                                 domain_index, &
                                 constraints,  &
                                 domain,       &
                                 enefunc)

    ! formal arguments
    type(s_par),          intent(in)    :: par
    type(s_prmtop),       intent(in)    :: prmtop
    type(s_grotop),       intent(in)    :: grotop
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_domain_index), intent(in)    :: domain_index
    type(s_constraints),  intent(inout) :: constraints
    type(s_domain),       intent(inout) :: domain
    type(s_enefunc),      intent(inout) :: enefunc


    ! charmm
    !
    if (par%num_bonds > 0) then

      call define_enefunc_charmm_peer (par, ene_info, domain_index, &
                                       constraints, domain, enefunc)

    ! amber
    !
    else if (prmtop%num_atoms > 0) then

      call define_enefunc_amber_peer  (prmtop, ene_info, domain_index, &
                                       constraints, domain, enefunc)

    ! gromacs
    !
    else if (grotop%num_atomtypes > 0) then

      call define_enefunc_gromacs_peer(grotop, ene_info, domain_index, &
                                       constraints, domain, enefunc)

    end if

    return

  end subroutine define_enefunc_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_charmm_peer
  !> @brief        peer routine of
  !!                            sp_enefunc_charmm_mod::define_enefunc_charmm()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_charmm_peer(par,          &
                                        ene_info,     &
                                        domain_index, &
                                        constraints,  &
                                        domain,       &
                                        enefunc)

    ! formal arguments
    type(s_par),          intent(in)    :: par
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_domain_index), intent(in)    :: domain_index
    type(s_constraints),  intent(inout) :: constraints
    type(s_domain),       intent(inout) :: domain
    type(s_enefunc),      intent(inout) :: enefunc

    ! local variables
    integer               :: ncel, ncelb



    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)


    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond_C_peer(par, domain_index, domain, enefunc)


      ! angle
      !
      call setup_enefunc_angl_C_peer(par, domain_index, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint_C_peer( &
                            par, constraints, domain_index, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint_C_peer( &
                            par, constraints, domain_index, domain, enefunc)

    end if


    ! dihedral
    !
    call setup_enefunc_dihe_C_peer(par, domain_index, domain, enefunc)


    ! improper
    !
    call setup_enefunc_impr_C_peer(par, domain_index, domain, enefunc)


    ! cmap
    !
    call setup_enefunc_cmap_C_peer(par, domain_index, domain, enefunc, &
                                        ene_info)

    ! nonbonded
    !
    call setup_enefunc_nonb_C_peer(par, constraints, domain, enefunc, ene_info)


    ! restraint
    !
    call setup_enefunc_restraints_peer(domain_index, domain, enefunc)

    return

  end subroutine define_enefunc_charmm_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_bond()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_C_peer(par, domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variable
    integer                       :: i, j, icel_local
    integer                       :: icel1, icel2
    integer                       :: nbond_p, i1, i2
    integer                       :: nbond_lst, ibond
    character(6)                  :: ci1, ci2

    real(wp),             pointer :: force(:,:), dist(:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: bond(:), list(:,:,:)
    integer,              pointer :: ibond_lst(:)


    ibond_lst => domain_index%ibond_list

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min

    nbond_lst = domain_index%nbond_list
    nbond_p   = par%num_bonds


    do ibond = 1, nbond_lst

      i = ibond_lst(ibond)

      i1  = hm_bond_list(1, i)
      i2  = hm_bond_list(2, i)
      ci1 = hm_atom_cls_name(i1)
      ci2 = hm_atom_cls_name(i2)

      icel1 = domain%id_g2l(1, i1)
      icel2 = domain%id_g2l(1, i2)

      ! Check if it is in my domain
      !
      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1, icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          do j = 1, nbond_p
            if ((ci1 == par%bond_atom_cls(1, j) .and.  &
                 ci2 == par%bond_atom_cls(2, j)) .or.  &
                (ci1 == par%bond_atom_cls(2, j) .and.  &
                 ci2 == par%bond_atom_cls(1, j))) then
              bond (icel_local) = bond(icel_local) + 1
              list (1, bond(icel_local), icel_local) = i1
              list (2, bond(icel_local), icel_local) = i2
              force(bond(icel_local), icel_local) = par%bond_force_const(j) 
              dist (bond(icel_local), icel_local) = par%bond_dist_min(j)
              exit
            end if
          end do

          if (j == nbond_p + 1) &
            write(MsgOut,*) &
                 'Setup_Enefunc_Bond_Peer> not found BOND: [', &
                 ci1, ']-[', ci2, '] in parameter file. (ERROR)'

        end if
      end if

    end do

    return

  end subroutine setup_enefunc_bond_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_angl()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_C_peer(par, domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                       :: i, j, icel_local
    integer                       :: icel1, icel2
    integer                       :: nangl_p
    integer                       :: list(3)
    integer                       :: nangl_lst, iangl
    character(6)                  :: ci1, ci2, ci3

    real(wp),             pointer :: force(:,:), theta(:,:)
    real(wp),             pointer :: ubforce(:,:), ubrmin(:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: angle(:), alist(:,:,:)
    integer,              pointer :: iangl_list(:)


    iangl_list => domain_index%iangl_list

    ncel       => domain%num_cell_local
    cell_pair  => domain%cell_pair

    angle      => enefunc%num_angle
    alist      => enefunc%angle_list
    force      => enefunc%angle_force_const
    theta      => enefunc%angle_theta_min
    ubforce    => enefunc%urey_force_const
    ubrmin     => enefunc%urey_rmin

    nangl_lst  = domain_index%nangl_list
    nangl_p    = par%num_angles


    do iangl = 1, nangl_lst

      i = iangl_list(iangl)

      list(1) = hm_angl_list(1, i)
      list(2) = hm_angl_list(2, i)
      list(3) = hm_angl_list(3, i)
      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))

      icel1 = domain%id_g2l(1, list(1))
      icel2 = domain%id_g2l(1, list(3))

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1, icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          do j = 1, nangl_p
            if ((ci1 == par%angl_atom_cls(1, j) .and. &
                 ci2 == par%angl_atom_cls(2, j) .and. &
                 ci3 == par%angl_atom_cls(3, j)) .or. &
                (ci1 == par%angl_atom_cls(3, j) .and. &
                 ci2 == par%angl_atom_cls(2, j) .and. &
                 ci3 == par%angl_atom_cls(1, j))) then
              angle(icel_local) = angle(icel_local) + 1
              alist(1:3,angle(icel_local),icel_local) = list(1:3)
              force  (angle(icel_local),icel_local) = par%angl_force_const(j)
              theta  (angle(icel_local),icel_local) = par%angl_theta_min(j)*RAD
              ubforce(angle(icel_local),icel_local) = par%angl_ub_force_const(j)
              ubrmin (angle(icel_local),icel_local) = par%angl_ub_rmin(j)
              exit
            end if
          end do

          if (j == nangl_p + 1) &
            write(MsgOut,*) &
              'Setup_Enefunc_Angl_Peer> not found ANGL: [', &
              ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'

        end if
      end if

    end do

    return

  end subroutine setup_enefunc_angl_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_dihe()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_C_peer(par, domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                       :: ndihe_p
    integer                       :: i, j, icel_local
    integer                       :: icel1, icel2
    integer                       :: nw_found
    integer                       :: list(4)
    integer                       :: ndihe_lst, idihe
    character(6)                  :: ci1, ci2, ci3, ci4

    real(wp),             pointer :: force(:,:), phase(:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: dihedral(:), dlist(:,:,:), period(:,:)
    integer,              pointer :: idihe_lst(:)
    integer,              pointer :: notation

    logical,          allocatable :: no_wild(:)


    idihe_lst  => domain_index%idihe_list

    ncel       => domain%num_cell_local
    cell_pair  => domain%cell_pair

    dihedral   => enefunc%num_dihedral
    dlist      => enefunc%dihe_list
    force      => enefunc%dihe_force_const
    period     => enefunc%dihe_periodicity
    phase      => enefunc%dihe_phase
    notation   => enefunc%notation_14types
    notation   = 100

    ndihe_lst  = domain_index%ndihe_list
    ndihe_p    = par%num_dihedrals


    ! check usage of wild card
    !
    allocate(no_wild(ndihe_p))

    do i = 1, ndihe_p

      if ((par%dihe_atom_cls(1,i) /= WildCard) .and. &
          (par%dihe_atom_cls(4,i) /= WildCard)) then
        ! A-B-C-D type
        no_wild(i) = .true.
      else
        ! X-B-C-D type
        no_wild(i) = .false.
      end if

    end do

    ! find number of interactions
    !
    do idihe = 1, ndihe_lst

      i = idihe_lst(idihe)

      list(1) = hm_dihe_list(1, i)
      list(2) = hm_dihe_list(2, i)
      list(3) = hm_dihe_list(3, i)
      list(4) = hm_dihe_list(4, i)
      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))
      ci4 = hm_atom_cls_name(list(4))

      icel1 = domain%id_g2l(1, list(1))
      icel2 = domain%id_g2l(1, list(4))


      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          nw_found = 0
          do j = 1, ndihe_p
            if (no_wild(j)) then
              if (((ci1 == par%dihe_atom_cls(1, j)) .and. &
                   (ci2 == par%dihe_atom_cls(2, j)) .and. &
                   (ci3 == par%dihe_atom_cls(3, j)) .and. &
                   (ci4 == par%dihe_atom_cls(4, j))) .or. &
                  ((ci1 == par%dihe_atom_cls(4, j)) .and. &
                   (ci2 == par%dihe_atom_cls(3, j)) .and. &
                   (ci3 == par%dihe_atom_cls(2, j)) .and. &
                   (ci4 == par%dihe_atom_cls(1, j)))) then
                nw_found = nw_found + 1
                dihedral(icel_local) = dihedral(icel_local) + 1 
                dlist (1:4,dihedral(icel_local),icel_local) = list(1:4)
                force (dihedral(icel_local),icel_local) =par%dihe_force_const(j)
                period(dihedral(icel_local),icel_local) =par%dihe_periodicity(j)
                phase (dihedral(icel_local),icel_local) =par%dihe_phase(j) * RAD
                if (period(dihedral(icel_local),icel_local) >  &
                enefunc%notation_14types) &
                call error_msg('Setup_Enefunc_Dihe_C_Peer> Too many periodicity.')
              end if
            end if
          end do

          if (nw_found == 0) then
            do j = 1, ndihe_p
              if (.not.no_wild(j)) then
                if (((ci2 == par%dihe_atom_cls(2, j)) .and. &
                     (ci3 == par%dihe_atom_cls(3, j))) .or. &
                    ((ci2 == par%dihe_atom_cls(3, j)) .and. &
                     (ci3 == par%dihe_atom_cls(2, j)))) then
                  dihedral(icel_local) = dihedral(icel_local) + 1
                  dlist (1:4,dihedral(icel_local),icel_local) = list(1:4)
                  force (dihedral(icel_local),icel_local) = &
                       par%dihe_force_const(j)
                  period(dihedral(icel_local),icel_local) = &
                       par%dihe_periodicity(j)
                  phase (dihedral(icel_local),icel_local) = &
                       par%dihe_phase(j) * RAD
                  if (period(dihedral(icel_local),icel_local) >  &
                  enefunc%notation_14types) &
                  call error_msg('Setup_Enefunc_Dihe_C_Peer> Too many periodicity.')
                end if
              end if
            end do
          end if

        end if
      end if
    end do

    deallocate(no_wild)

    return

  end subroutine setup_enefunc_dihe_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_impr()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_C_peer(par, domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                       :: nimpr_p
    integer                       :: i, j, icel_local
    integer                       :: icel1, icel2
    integer                       :: list(4)
    integer                       :: nimpr_lst, iimpr
    character(6)                  :: ci1, ci2, ci3, ci4
    logical                       :: no_wild

    real(wp),             pointer :: force(:,:), phase(:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: improper(:), ilist(:,:,:)
    integer,              pointer :: iimpr_lst(:)
    integer,          allocatable :: wc_type(:)


    iimpr_lst  => domain_index%iimpr_list

    ncel       => domain%num_cell_local
    cell_pair  => domain%cell_pair

    improper   => enefunc%num_improper
    ilist      => enefunc%impr_list
    force      => enefunc%impr_force_const
    phase      => enefunc%impr_phase

    nimpr_lst  = domain_index%nimpr_list
    nimpr_p    = par%num_impropers


    ! check usage of wild card
    !
    allocate(wc_type(nimpr_p))

    do i = 1, nimpr_p

      if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
          (par%impr_atom_cls(2,i) /= WildCard) .and. &
          (par%impr_atom_cls(3,i) /= WildCard) .and. &
          (par%impr_atom_cls(4,i) /= WildCard)) then
        ! A-B-C-D type
        wc_type(i) = 0

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) /= WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! X-B-C-D type
        wc_type(i) = 1

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! X-X-C-D type
        wc_type(i) = 2

      else if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) == WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then
        ! A-X-X-D type
        wc_type(i) = 3

      else
        call error_msg('Setup_Enefunc_Impr_Peer> Undefined Wild Card')

      end if

    end do

    ! setup parameters
    !

    do iimpr = 1, nimpr_lst

      i = iimpr_lst(iimpr)

      no_wild = .false.

      list(1) = hm_impr_list(1, i)
      list(2) = hm_impr_list(2, i)
      list(3) = hm_impr_list(3, i)
      list(4) = hm_impr_list(4, i)
      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))
      ci4 = hm_atom_cls_name(list(4))

      icel1 = domain%id_g2l(1, list(1))
      icel2 = domain%id_g2l(1, list(4))

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          ! A-B-C-D type
          !
          do j = 1, nimpr_p
            if (wc_type(j) == 0) then
              if (((ci1 == par%impr_atom_cls(1, j)) .and. &
                   (ci2 == par%impr_atom_cls(2, j)) .and. &
                   (ci3 == par%impr_atom_cls(3, j)) .and. &
                   (ci4 == par%impr_atom_cls(4, j))) .or. &
                  ((ci1 == par%impr_atom_cls(4, j)) .and. &
                   (ci2 == par%impr_atom_cls(3, j)) .and. &
                   (ci3 == par%impr_atom_cls(2, j)) .and. &
                   (ci4 == par%impr_atom_cls(1, j)))) then

                improper(icel_local) = improper(icel_local) + 1
                ilist(1:4,improper(icel_local), icel_local) = list(1:4)
                force(improper(icel_local),icel_local) = par%impr_force_const(j)
                phase(improper(icel_local),icel_local) = par%impr_phase(j) * RAD
                no_wild = .true.
                exit

              end if
            end if
          end do

          ! X-B-C-D type
          !
          if (.not.no_wild) then
            do j = 1, nimpr_p
              if (wc_type(j) == 1) then
                if (((ci2 == par%impr_atom_cls(2, j)) .and. &
                     (ci3 == par%impr_atom_cls(3, j)) .and. &
                     (ci4 == par%impr_atom_cls(4, j))) .or. &
                    ((ci2 == par%impr_atom_cls(4, j)) .and. &
                     (ci3 == par%impr_atom_cls(3, j)) .and. &
                     (ci4 == par%impr_atom_cls(2, j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local), icel_local) = list(1:4)
                  force(improper(icel_local),icel_local)=par%impr_force_const(j)
                  phase(improper(icel_local),icel_local)=par%impr_phase(j)*RAD
                  no_wild = .true.
                  exit

                end if
              end if
            end do
          end if

          ! X-X-C-D type
          !
          if (.not.no_wild) then
            do j = 1, nimpr_p
              if (wc_type(j) == 2) then
                if (((ci3 == par%impr_atom_cls(3, j)) .and. &
                     (ci4 == par%impr_atom_cls(4, j))) .or. &
                    ((ci3 == par%impr_atom_cls(4, j)) .and. &
                     (ci4 == par%impr_atom_cls(3, j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local), icel_local) = list(1:4)
                  force(improper(icel_local),icel_local)=par%impr_force_const(j)
                  phase(improper(icel_local),icel_local)=par%impr_phase(j)*RAD
                  no_wild = .true.
                  exit

                end if
              end if
            end do
          end if

          ! A-X-X-D type
          !
          if (.not.no_wild) then
            do j = 1, nimpr_p
              if (wc_type(j) == 3) then
                if (((ci1 == par%impr_atom_cls(1, j)) .and. &
                     (ci4 == par%impr_atom_cls(4, j))) .or. &
                    ((ci1 == par%impr_atom_cls(4, j)) .and. &
                     (ci4 == par%impr_atom_cls(1, j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local), icel_local) = list(1:4)
                  force(improper(icel_local),icel_local)=par%impr_force_const(j)
                  phase(improper(icel_local),icel_local)=par%impr_phase(j)*RAD
                  no_wild = .true.
                  exit

                end if
              end if
            end do
          end if

          if (.not.no_wild) &
            write(MsgOut,*) &
              'Setup_Enefunc_Impr_Peer> Unknown IMPR type. [', &
              ci1, ']-[', ci2, ']-[', ci3, ']-[', ci4, '] (ERROR)'

        end if
      end if

    end do

    deallocate(wc_type)

    return

  end subroutine setup_enefunc_impr_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_cmap()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap_C_peer(par, domain_index, &
                                       domain, enefunc, ene_info)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc
    type(s_ene_info),             intent(in)    :: ene_info

    ! local variables
    integer                       :: ncmap_p, ngrid0
    integer                       :: list(2), icel1, icel2, icel_local
    integer                       :: flag_cmap_type, alloc_stat, dealloc_stat
    integer                       :: i, j, k, l, ityp
    integer                       :: ncmap_lst, icmap
    character(6)                  :: ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8
    logical                       :: periodic

    real(wp),             pointer :: ccoef(:,:,:,:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: ncmap(:), clist(:,:,:), ctype(:,:), cres(:)
    integer,              pointer :: icmap_lst(:)

    real(wp),         allocatable :: c_ij(:,:,:,:)      ! cmap coeffs


    icmap_lst  => domain_index%icmap_list

    ncel       => domain%num_cell_local
    cell_pair  => domain%cell_pair

    ! If 'periodic' is .true.,
    !   then cubic spline with periodic (in dihedral-angle space) boundary
    !   will be applied.
    ! If 'periodic' is .false.,
    !   then natural cubic spline without periodic boudnary
    !   will be applied to expanded cross-term maps.
    !   This is similar with CHARMM's source code.
    !
    periodic   = ene_info%cmap_pspline
    ncmap_lst  = domain_index%ncmap_list
    ncmap_p    = par%num_cmaps

    !  begin
    !
    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, par%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncel, ngrid0, ncmap_p)

    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ncmap      => enefunc%num_cmap
    clist      => enefunc%cmap_list
    ctype      => enefunc%cmap_type
    cres       => enefunc%cmap_resolution
    ccoef      => enefunc%cmap_coef

    do i = 1, ncmap_p
      cres(i) = par%cmap_resolution(i)
    end do

    !  derive cmap coefficients by bicubic interpolation
    !
    do ityp = 1, ncmap_p
      if (periodic) then
        call derive_cmap_coefficients_p(ityp, par, c_ij)
      else
        call derive_cmap_coefficients_np(ityp, par, c_ij)
      end if
      do l = 1, ngrid0
        do k = 1, ngrid0 
          do j = 1, 4
            do i = 1, 4
              ccoef(i,j,k,l,ityp) = c_ij(i,j,k,l)
            end do
          end do
        end do
      end do
    end do

    ncmap(1:ncel) = 0

    do icmap = 1, ncmap_lst

      i = icmap_lst(icmap)

      list(1) = hm_cmap_list(1, i)
      list(2) = hm_cmap_list(8, i)
      icel1 = domain%id_g2l(1, list(1))
      icel2 = domain%id_g2l(1, list(2))

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = domain%cell_pair(icel1, icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          ! ci* will be atom-type strings
          !
          ci1 = hm_atom_cls_name(hm_cmap_list(1,i))
          ci2 = hm_atom_cls_name(hm_cmap_list(2,i))
          ci3 = hm_atom_cls_name(hm_cmap_list(3,i))
          ci4 = hm_atom_cls_name(hm_cmap_list(4,i))
          ci5 = hm_atom_cls_name(hm_cmap_list(5,i))
          ci6 = hm_atom_cls_name(hm_cmap_list(6,i))
          ci7 = hm_atom_cls_name(hm_cmap_list(7,i))
          ci8 = hm_atom_cls_name(hm_cmap_list(8,i))
          flag_cmap_type = -1

          ! assign cmap type ID to each (psi,phi) pair
          !
          do j = 1, ncmap_p
            if (ci1 == par%cmap_atom_cls(1, j) .and. &
                ci2 == par%cmap_atom_cls(2, j) .and. &
                ci3 == par%cmap_atom_cls(3, j) .and. &
                ci4 == par%cmap_atom_cls(4, j) .and. &
                ci5 == par%cmap_atom_cls(5, j) .and. &
                ci6 == par%cmap_atom_cls(6, j) .and. &
                ci7 == par%cmap_atom_cls(7, j) .and. &
                ci8 == par%cmap_atom_cls(8, j)) then

              ncmap(icel_local) = ncmap(icel_local) + 1
              clist(1, ncmap(icel_local), icel_local) = hm_cmap_list(1, i)
              clist(2, ncmap(icel_local), icel_local) = hm_cmap_list(2, i)
              clist(3, ncmap(icel_local), icel_local) = hm_cmap_list(3, i)
              clist(4, ncmap(icel_local), icel_local) = hm_cmap_list(4, i)
              clist(5, ncmap(icel_local), icel_local) = hm_cmap_list(5, i)
              clist(6, ncmap(icel_local), icel_local) = hm_cmap_list(6, i)
              clist(7, ncmap(icel_local), icel_local) = hm_cmap_list(7, i)
              clist(8, ncmap(icel_local), icel_local) = hm_cmap_list(8, i)
              ctype(ncmap(icel_local), icel_local) = j
              flag_cmap_type = j
              exit
            end if
          end do

          ! if not found, print detail about the error.
          !
          if (flag_cmap_type <= 0) then
            write(MsgOut,*) 'Setup_Enefunc_Cmap_Peer> not found CMAP: '
            write(MsgOut,*) ' [',ci1,']-[',ci2,']-[',ci3,']-[',ci4,']-'
            write(MsgOut,*) ' [',ci5,']-[',ci6,']-[',ci7,']-[',ci8,'] '
            write(MsgOut,*)  'in parameter file. (ERROR)'
          end if

        end if
      end if

    end do

    deallocate(c_ij)

    return

  end subroutine setup_enefunc_cmap_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::
  !!                                          setup_enefunc_bond_constraint()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint_C_peer(par, constraints, &
                                                  domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_constraints),  target, intent(inout) :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variable
    integer                       :: i, j, k, m, ih, icel_local, connect
    integer                       :: icel1, icel2, icel
    integer                       :: nbond_p, i1, i2, ih1, ih2
    integer                       :: nbond_lst, ibond, nbond_c
    character(6)                  :: ci1, ci2
    character(4)                  :: a1, a2

    real(wp),             pointer :: force(:,:), dist(:,:)
    real(wp),             pointer :: HGr_bond_dist(:,:,:,:)
    integer,              pointer :: bond(:), list(:,:,:), ncel
    integer,              pointer :: cell_pair(:,:), id_g2l(:,:), id_l2g(:,:)
    integer,              pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,              pointer :: ibond_lst(:)


    ibond_lst     => domain_index%ibond_list

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min

    nbond_lst     =  domain_index%nbond_list
    nbond_p       =  par%num_bonds

    connect       = constraints%connect

    bond(1:ncel)  = 0

    nbond_c       = 0

    do ibond = 1, nbond_lst

      i = ibond_lst(ibond)

      i1 = hm_bond_list(1,i)
      i2 = hm_bond_list(2,i)

      ci1 = hm_atom_cls_name(i1)
      ci2 = hm_atom_cls_name(i2)

      a1 = hm_atom_name(i1)
      a2 = hm_atom_name(i2)

      if (a1(1:1) /= 'H' .and. a2(1:1) /= 'H' .and. &
          a1(1:1) /= 'h' .and. a2(1:1) /= 'h') then

        icel1 = domain%id_g2l(1,i1)
        icel2 = domain%id_g2l(1,i2)

        ! Check if it is in my domain
        !
        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1, icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            do j = 1, nbond_p
              if ((ci1 == par%bond_atom_cls(1, j) .and.  &
                   ci2 == par%bond_atom_cls(2, j)) .or.  &
                  (ci1 == par%bond_atom_cls(2, j) .and.  &
                   ci2 == par%bond_atom_cls(1, j))) then
                bond (icel_local) = bond(icel_local) + 1
                list (1, bond(icel_local), icel_local) = i1
                list (2, bond(icel_local), icel_local) = i2
                force(bond(icel_local), icel_local) = par%bond_force_const(j) 
                dist (bond(icel_local), icel_local) = par%bond_dist_min(j)
                exit
              end if
            end do

            if (j == nbond_p + 1) &
              write(MsgOut,*) &
                'Setup_Enefunc_Bond_Constraint_Peer> not found BOND: [', &
                ci1, ']-[', ci2, '] in parameter file. (ERROR)'

          end if

        end if

      else

        do icel = 1, ncel
          do j = 1, connect

            do k = 1, HGr_local(j,icel)
              ih1 = id_l2g(HGr_bond_list(1,k,j,icel),icel)
              do ih = 1, j
                ih2 = id_l2g(HGr_bond_list(ih+1,k,j,icel),icel)

                if (ih1 == i1 .and. ih2 == i2 .or. &
                    ih2 == i1 .and. ih1 == i2) then

                  do m = 1, nbond_p
                    if ((ci1 == par%bond_atom_cls(1, m) .and.  &
                         ci2 == par%bond_atom_cls(2, m)) .or.  &
                        (ci1 == par%bond_atom_cls(2, m) .and.  &
                         ci2 == par%bond_atom_cls(1, m))) then

                      nbond_c = nbond_c + 1
                      HGr_bond_dist(ih+1,k,j,icel) = par%bond_dist_min(m)
                      exit
                    end if
                  end do

                end if

              end do
            end do

          end do
        end do

      end if

    end do
    constraints%num_bonds = nbond_c

    return

  end subroutine setup_enefunc_bond_constraint_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::
  !!                                          setup_enefunc_angl_constraint()
  !! @authors      NT, JJ
  !! @param[in]    
  !! @date         2012/10/24 (JJ) - removing angle list of water molecules
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint_C_peer(par, constraints, &
                                                  domain_index, domain, enefunc)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_constraints),          intent(in)    :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                       :: i, j, icel_local
    integer                       :: icel1, icel2
    integer                       :: nangl_p
    integer                       :: list(3)
    integer                       :: nangl_lst, iangl
    character(6)                  :: ci1, ci2, ci3
    character(6)                  :: ri1, ri2, ri3

    real(wp),             pointer :: force(:,:), theta(:,:)
    real(wp),             pointer :: ubforce(:,:), ubrmin(:,:)
    integer,              pointer :: ncel, cell_pair(:,:)
    integer,              pointer :: angle(:), alist(:,:,:)
    integer,              pointer :: iangl_list(:)


    iangl_list => domain_index%iangl_list

    ncel       => domain%num_cell_local
    cell_pair  => domain%cell_pair

    angle      => enefunc%num_angle
    alist      => enefunc%angle_list
    force      => enefunc%angle_force_const
    theta      => enefunc%angle_theta_min
    ubforce    => enefunc%urey_force_const
    ubrmin     => enefunc%urey_rmin

    nangl_lst  = domain_index%nangl_list
    nangl_p    = par%num_angles


    do iangl = 1, nangl_lst

      i = iangl_list(iangl)

      list(1) = hm_angl_list(1, i)
      list(2) = hm_angl_list(2, i)
      list(3) = hm_angl_list(3, i)
      ci1 = hm_atom_cls_name(list(1))
      ci2 = hm_atom_cls_name(list(2))
      ci3 = hm_atom_cls_name(list(3))
      ri1 = hm_residue_name(list(1))
      ri2 = hm_residue_name(list(2))
      ri3 = hm_residue_name(list(3))

      if (ri1(1:4) /= constraints%water_model .and. &
          ri2(1:4) /= constraints%water_model .and. &
          ri3(1:4) /= constraints%water_model) then

        icel1 = domain%id_g2l(1, list(1))
        icel2 = domain%id_g2l(1, list(3))

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1, icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            do j = 1, nangl_p
              if ((ci1 == par%angl_atom_cls(1, j) .and. &
                   ci2 == par%angl_atom_cls(2, j) .and. &
                   ci3 == par%angl_atom_cls(3, j)) .or. &
                  (ci1 == par%angl_atom_cls(3, j) .and. &
                   ci2 == par%angl_atom_cls(2, j) .and. &
                   ci3 == par%angl_atom_cls(1, j))) then

                angle(icel_local) = angle(icel_local) + 1
                alist(1:3,angle(icel_local),icel_local) = list(1:3)
              force  (angle(icel_local),icel_local) = par%angl_force_const(j)
              theta  (angle(icel_local),icel_local) = par%angl_theta_min(j)*RAD
              ubforce(angle(icel_local),icel_local) = par%angl_ub_force_const(j)
              ubrmin (angle(icel_local),icel_local) = par%angl_ub_rmin(j)
                exit

              end if
            end do

            if (j == nangl_p + 1) &
              write(MsgOut,*) &
                'Setup_Enefunc_Angl_Constraint_Peer> not found ANGL: [', &
                ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'

          end if

        end if

      end if

    end do

    return

  end subroutine setup_enefunc_angl_constraint_C_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_C_peer
  !> @brief        peer routine of sp_enefunc_charmm_mod::setup_enefunc_nonb()
  !! @authors      NT (modified by JJ)
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_C_peer(par, constraints, domain,  &
                                       enefunc, ene_info)

    ! formal arguments
    type(s_par),                  intent(in)    :: par
    type(s_constraints),          intent(in)    :: constraints
    type(s_domain),       target, intent(inout) :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc
    type(s_ene_info),             intent(in)    :: ene_info

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin
    integer                  :: i, j, k, ix, jx, kx
    integer                  :: identical_sign
    integer                  :: nonb_p, nbfx_p, cls_local, ncel
    character(6)             :: ci1, ci2

    integer,  allocatable    :: nonb_atom_cls(:), check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)

    enefunc%num_atom_cls = par%num_atom_cls

    ELECOEF          = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    nonb_p = enefunc%num_atom_cls

    allocate(nonb_atom_cls(nonb_p),         &
             check_cls(nonb_p),             &
             atmcls_map_g2l(nonb_p),        &
             atmcls_map_l2g(nonb_p),        &
             nb14_lj6     (nonb_p, nonb_p), &
             nb14_lj12    (nonb_p, nonb_p), &
             nonb_lj6     (nonb_p, nonb_p), &
             nonb_lj12    (nonb_p, nonb_p))

    nonb_atom_cls(1:nonb_p)           = 0
    check_cls(1:nonb_p)               = 0
    nb14_lj6     (1:nonb_p, 1:nonb_p) = 0.0_wp
    nb14_lj12    (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj6     (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj12    (1:nonb_p, 1:nonb_p) = 0.0_wp

    do i = 1, nonb_p

      lj_coef(1,i) = par%nonb_eps(i)
      lj_coef(2,i) = par%nonb_rmin(i)

      do j = 1, nonb_p

        ! combination rule
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
        nb14_lj6(i,j)  = 2.0_wp * eps14 * (rmin14 ** 6)
        nonb_lj12(i,j) = eps * (rmin ** 12)
        nonb_lj6(i,j)  = 2.0_wp * eps * (rmin ** 6)

      end do
    end do

    ! overwrite lennard-jones parameters by NBFIX parameters
    !
    nbfx_p = par%num_nbfix

    do k = 1, nbfx_p
      ci1 = par%nbfi_atom_cls(1,k)
      ci2 = par%nbfi_atom_cls(2,k)
      do i = 1, nonb_p
        do j = 1, nonb_p
          if ((ci1 == par%nonb_atom_cls(i)  .and. &
               ci2 == par%nonb_atom_cls(j)) .or.  &
              (ci2 == par%nonb_atom_cls(i) .and. &
               ci1 == par%nonb_atom_cls(j))) then

            ! combination rule
            !
            eps14  = abs(par%nbfi_eps_14 (k)) !TODO CHECK
            rmin14 = par%nbfi_rmin_14(k)
            eps    = abs(par%nbfi_eps    (k))
            rmin   = par%nbfi_rmin   (k)

            ! set parameters
            !
            nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
            nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14 ** 6)
            nonb_lj12(i,j) = eps * (rmin ** 12)
            nonb_lj6 (i,j) = 2.0_wp * eps * (rmin ** 6)
          end if
        end do
      end do
    end do

    ! check the usage of atom class
    !
    do i = 1, hm_num_atoms
      k = hm_atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(hm_atom_cls_name(i))//'"')
      endif

      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nonb_p
      if (check_cls(i) == 1) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do

    domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
    domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
    domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
    domain%water%atom_cls_no(1:3) = &
      atmcls_map_g2l(domain%water%atom_cls_no(1:3))
    enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
    enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
    if (constraints%tip4 .or. enefunc%table%tip4) then
      enefunc%table%atom_cls_no_D = atmcls_map_g2l(enefunc%table%atom_cls_no_D)
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
    end if

    deallocate(nonb_atom_cls, &
               check_cls,     &
               atmcls_map_g2l,&
               atmcls_map_l2g,&
               nb14_lj6,      &
               nb14_lj12,     &
               nonb_lj6,      &
               nonb_lj12)

    enefunc%num_atom_cls = cls_local

    return

  end subroutine setup_enefunc_nonb_C_peer


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber_peer
  !> @brief        peer routine of
  !!                            sp_enefunc_amber_mod::define_enefunc_amber()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber_peer(prmtop,       &
                                       ene_info,     &
                                       domain_index, &
                                       constraints,  &
                                       domain,       &
                                       enefunc)

    ! formal arguments
    type(s_prmtop),       intent(in)    :: prmtop
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_domain_index), intent(in)    :: domain_index
    type(s_constraints),  intent(inout) :: constraints
    type(s_domain),       intent(inout) :: domain
    type(s_enefunc),      intent(inout) :: enefunc

    ! local variables
    integer               :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)


    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond_A_peer(prmtop, domain_index, domain, enefunc)


      ! angle
      !
      call setup_enefunc_angl_A_peer(prmtop, domain_index, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint_A_peer( &
                           prmtop, constraints, domain_index, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint_A_peer( &
                           prmtop, constraints, domain_index, domain, enefunc)

    end if


    ! dihedral
    !
    call setup_enefunc_dihe_A_peer(prmtop, domain_index, domain, enefunc)


    ! improper
    !
    call setup_enefunc_impr_A_peer(prmtop, domain_index, domain, enefunc)


    ! nonbonded
    !
    call setup_enefunc_nonb_A_peer(ene_info, prmtop, constraints,  &
                                   domain, enefunc) 


    ! restraint
    !
    call setup_enefunc_restraints_peer(domain_index, domain, enefunc)

    return

  end subroutine define_enefunc_amber_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_bond()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_A_peer(prmtop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, i1, i2, icel1, icel2
    integer                  :: icel_local

    real(wp),        pointer :: force(:,:), dist(:,:)
    integer,         pointer :: bond(:), list(:,:,:)
    integer,         pointer :: ncel, cell_pair(:,:)
    integer,         pointer :: id_g2l(:,:)
    integer,         pointer :: nbond


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min

    do i = 1, prmtop%num_bondh

      i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
      i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i2)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nbond => bond(icel_local)
          nbond = nbond + 1

          if (nbond > MaxBond) &
            call error_msg('Setup_Enefunc_Bond_A_Peer> Too many bonds.')

          list (1,nbond,icel_local) = i1
          list (2,nbond,icel_local) = i2
          force(  nbond,icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
          dist (  nbond,icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
        end if

      end if

    end do

    do i = 1, prmtop%num_mbonda

      i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
      i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i2)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nbond => bond(icel_local)
          nbond = nbond + 1

          if (nbond > MaxBond) &
            call error_msg('Setup_Enefunc_Bond_A_Peer> Too many bonds.')

          list (1,nbond,icel_local) = i1
          list (2,nbond,icel_local) = i2
          force(  nbond,icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
          dist (  nbond,icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
        end if

      end if

    end do

    return

  end subroutine setup_enefunc_bond_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_angl()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_A_peer(prmtop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, i1, i2, i3, icel1, icel2
    integer                  :: icel_local

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    integer,         pointer :: nangl


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min


    do i = 1, prmtop%num_anglh

      i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
      i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
      i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i3)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nangl => angle(icel_local)
          nangl = nangl + 1

          if (nangl > MaxAngle) &
            call error_msg('Setup_Enefunc_Angl_A_Peer> Too many angles.')

          alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
          force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
          theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
        end if

      end if

    end do

    do i = 1, prmtop%num_mangla

      i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
      i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
      i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i3)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nangl => angle(icel_local)
          nangl = nangl + 1

          if (nangl > MaxAngle) &
            call error_msg('Setup_Enefunc_Angl_A_Peer> Too many angles.')

          alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
          force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
          theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
        end if

      end if

    end do

    return

  end subroutine setup_enefunc_angl_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_dihe()
  !! @authors      NT
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_A_peer(prmtop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                     :: i
    integer                     :: i1, i2, i3, i4, icel1, icel2, icel_local

    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types

    notation = 100
    if (prmtop%num_uniqdihe > 100) then
      notation = 1000
      if (prmtop%num_uniqdihe > 1000) then
        call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 
      endif
    endif

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
      i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
      i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
      i4 =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i4)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          ndihe => dihe(icel_local)
          ndihe = ndihe + 1

          if (ndihe > MaxDihe) &
            call error_msg('Setup_Enefunc_Dihe_A_Peer> Too many dihedrals.')

          list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
          force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
          period(  ndihe,icel_local) = &
                              prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
          if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
            period(ndihe,icel_local) = period(  ndihe,icel_local) &
                            + prmtop%dihe_inc_hy(5,i)*notation
          endif

        end if

      end if

    end do

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
      i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
      i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
      i4 =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i4)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          ndihe => dihe(icel_local)
          ndihe = ndihe + 1

          if (ndihe > MaxDihe) &
            call error_msg('Setup_Enefunc_Dihe_A_Peer> Too many dihedrals.')

          list(1:4,ndihe,icel_local) = (/i1,i2,i3,i4/)
          force (  ndihe,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          phase (  ndihe,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
          period(  ndihe,icel_local) = &
                              prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
          if (prmtop%lscee_scale_factor .or. prmtop%lscnb_scale_factor ) then
            period(ndihe,icel_local) = period(  ndihe,icel_local) &
                            + prmtop%dihe_wo_hy(5,i)*notation
          endif

        end if

      end if

    end do

    return

  end subroutine setup_enefunc_dihe_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_impr()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_A_peer(prmtop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                     :: i
    integer                     :: i1, i2, i3, i4, icel1, icel2, icel_local

    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) >= 0) &
        cycle

      i1 =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
      i2 =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
      i3 = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
      i4 = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i4)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nimpr => impr(icel_local)
          nimpr = nimpr + 1

          if (nimpr > MaxImpr) &
            call error_msg('Setup_Enefunc_Impr_A_Peer> Too many impropers.')

          list(1:4,nimpr,icel_local) = (/i1,i2,i3,i4/)
          force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
          phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))
          period(  nimpr,icel_local) = &
                              prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))

        end if

      end if

    end do

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) >= 0) &
        cycle

      i1 =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
      i2 =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
      i3 = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
      i4 = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i4)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nimpr => impr(icel_local)
          nimpr = nimpr + 1

          if (nimpr > MaxImpr) &
            call error_msg('Setup_Enefunc_Impr_A_Peer> Too many impropers.')

          list(1:4,nimpr,icel_local) = (/i1,i2,i3,i4/)
          force (  nimpr,icel_local) = &
                              prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
          phase (  nimpr,icel_local) = &
                              prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))
          period(  nimpr,icel_local) = &
                              prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))

        end if

      end if

    end do

    return

  end subroutine setup_enefunc_impr_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::
  !!                                          setup_enefunc_bond_constraint()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint_A_peer(prmtop, constraints, &
                                                 domain_index, domain,  enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_constraints),  target, intent(inout) :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                      :: i, j, k, m, ih, icel_local, connect
    integer                      :: i1, i2, ih1, ih2, icel1, icel2, icel
    integer                      :: wat_bonds, wat_found, nbond_c
    integer                      :: tmp_mole_no, mole_no
    character(6)                 :: ri1, ri2

    real(wp),            pointer :: force(:,:), dist(:,:)
    real(wp),            pointer :: HGr_bond_dist(:,:,:,:)
    integer,             pointer :: bond(:), list(:,:,:), ncel
    integer,             pointer :: cell_pair(:,:), id_g2l(:,:), id_l2g(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)


    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min

    connect       =  constraints%connect
    nbond_c       = 0

    do i = 1, prmtop%num_mbonda

      i1 = prmtop%bond_wo_hy(1,i) / 3 + 1
      i2 = prmtop%bond_wo_hy(2,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i2)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          bond(icel_local) = bond(icel_local) + 1

          if (bond(icel_local) > MaxBond) &
       call error_msg('Setup_Enefunc_Bond_Constraint_A_Peer> Too many bonds.')

          list (1,bond(icel_local),icel_local) = i1
          list (2,bond(icel_local),icel_local) = i2
          force(  bond(icel_local),icel_local) = &
                             prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
          dist (  bond(icel_local),icel_local) = &
                             prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
        end if

      end if

    end do

    do i = 1, prmtop%num_bondh

      i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
      i2 = prmtop%bond_inc_hy(2,i) / 3 + 1

      do icel = 1, ncel
        do j = 1, connect

          do k = 1, HGr_local(j,icel)
            ih1 = id_l2g(HGr_bond_list(1,k,j,icel),icel)
            do ih = 1, j
              ih2 = id_l2g(HGr_bond_list(ih+1,k,j,icel),icel)

              if (ih1 == i1 .and. ih2 == i2 .or. &
                  ih2 == i1 .and. ih1 == i2) then

                nbond_c = nbond_c + 1
                HGr_bond_dist(ih+1,k,j,icel) = &
                             prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))

              end if

            end do
          end do
        end do
      end do

    end do

    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      tmp_mole_no = -1

      do i = 1, prmtop%num_bondh
  
        i1 = prmtop%bond_inc_hy(1,i) / 3 + 1
        i2 = prmtop%bond_inc_hy(2,i) / 3 + 1
  
        ri1 = hm_residue_name(i1)
        ri2 = hm_residue_name(i2)

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds = wat_bonds+1
          mole_no = hm_molecule_no(hm_bond_list(1,i))

          if (mole_no /= tmp_mole_no) then
            wat_found = wat_found +1
            tmp_mole_no = mole_no
          end if

        end if

      end do

      if (wat_found /= enefunc%table%num_water) &
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint_A_Peer> # of water is incorrect')

    end if
    constraints%num_bonds = nbond_c

    return

  end subroutine setup_enefunc_bond_constraint_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::
  !!                                          setup_enefunc_angl_constraint()
  !! @authors      NT, JJ
  !! @param[in]    
  !! @date         2012/10/24 (JJ) - removing angle list of water molecules
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint_A_peer(prmtop, constraints, &
                                                  domain_index, domain, enefunc)

    ! formal arguments
    type(s_prmtop),               intent(in)    :: prmtop
    type(s_constraints),          intent(in)    :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, i1, i2, i3, icel1, icel2
    integer                  :: icel_local
    character(6)             :: res1, res2, res3, res4

    real(wp),        pointer :: force(:,:), theta(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    integer,         pointer :: nangl


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min


    do i = 1, prmtop%num_anglh

      i1 = prmtop%angl_inc_hy(1,i) / 3 + 1
      i2 = prmtop%angl_inc_hy(2,i) / 3 + 1
      i3 = prmtop%angl_inc_hy(3,i) / 3 + 1

      res1 = hm_residue_name(i1)
      res2 = hm_residue_name(i2)
      res3 = hm_residue_name(i3)

      if (res1 /= constraints%water_model .and. &
          res2 /= constraints%water_model .and. &
          res3 /= constraints%water_model) then

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i3)
  
        if (icel1 /= 0 .and. icel2 /= 0) then
  
          icel_local = cell_pair(icel1,icel2)
  
          if (icel_local > 0 .and. icel_local <= ncel) then
  
            nangl => angle(icel_local)
            nangl = nangl + 1
  
            if (nangl > MaxAngle) &
       call error_msg('Setup_Enefunc_Angl_Constraint_A_Peer> Too many angles.') 
  
            alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
            force(    nangl,icel_local) = &
                               prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
            theta(    nangl,icel_local) = &
                               prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))
          end if
  
        end if

      end if

    end do

    do i = 1, prmtop%num_mangla

      i1 = prmtop%angl_wo_hy(1,i) / 3 + 1
      i2 = prmtop%angl_wo_hy(2,i) / 3 + 1
      i3 = prmtop%angl_wo_hy(3,i) / 3 + 1

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i3)

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          nangl => angle(icel_local)
          nangl = nangl + 1

          if (nangl > MaxAngle) &
      call error_msg('Setup_Enefunc_Angl_Constraint_A_Peer> Too many angles.') 

          alist(1:3,nangl,icel_local) = (/i1, i2, i3/)
          force(    nangl,icel_local) = &
                             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
          theta(    nangl,icel_local) = &
                             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))
        end if

      end if

    end do

    return

  end subroutine setup_enefunc_angl_constraint_A_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_A_peer
  !> @brief        peer routine of sp_enefunc_amber_mod::setup_enefunc_nonb()
  !! @authors      NT (modified by JJ)
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_A_peer(ene_info, prmtop, constraints,  &
                                       domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, nnonb, ncel, ndihetype, ndihe
    integer                  :: nb_par_idx
    integer                  :: k, ix, jx, kx, cls_local

    integer,  allocatable    :: nonb_atom_cls(:), check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)

    ndihetype = size(prmtop%scee_scale_fact)
    call alloc_enefunc(enefunc, EneFuncAMBERScale, ndihetype, 0)

    ELECOEF          = ELECOEF_AMBER

    ndihe = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_inc_hy(5,i)) = &
         1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      ndihe = ndihe + 1

      if (prmtop%lscee_scale_factor) then
        enefunc%dihe_scee(prmtop%dihe_wo_hy(5,i)) =  &
          1.0_wp/prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scee(0) = 1.0_wp/1.2_wp
      end if

      if (prmtop%lscnb_scale_factor) then
        enefunc%dihe_scnb(prmtop%dihe_wo_hy(5,i)) = &
          1.0_wp/prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
      else
        enefunc%dihe_scnb(0) = 1.0_wp/2.0_wp
      end if

    end do

    nnonb = prmtop%num_types

    allocate(check_cls(nnonb),             &
             atmcls_map_g2l(nnonb),        &
             atmcls_map_l2g(nnonb),        &
             nb14_lj6      (nnonb, nnonb), &
             nb14_lj12     (nnonb, nnonb), &
             nonb_lj6      (nnonb, nnonb), &
             nonb_lj12     (nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          nb14_lj12(i,j) = 0.0_wp
          nb14_lj6 (i,j) = 0.0_wp
          nonb_lj12(i,j) = 0.0_wp
          nonb_lj6 (i,j) = 0.0_wp
        else
          nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

      end do
    end do

    ! check the usage of atom class
    !
    do i = 1, hm_num_atoms
      k = hm_atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(hm_atom_cls_name(i))//'"')
      endif
      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nnonb
      if (check_cls(i) == 1) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do
    if (constraints%rigid_bond .or. enefunc%table%water_table) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(1:3) = &
        atmcls_map_g2l(domain%water%atom_cls_no(1:3))
      enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
    endif
    if (constraints%tip4 .or. enefunc%table%tip4) then
      enefunc%table%atom_cls_no_D = atmcls_map_g2l(enefunc%table%atom_cls_no_D)
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
    end if

    deallocate(check_cls,     &
               atmcls_map_g2l,&
               atmcls_map_l2g,&
               nb14_lj6,      &
               nb14_lj12,     &
               nonb_lj6,      &
               nonb_lj12)

    enefunc%num_atom_cls = cls_local

    return

  end subroutine setup_enefunc_nonb_A_peer


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs_peer
  !> @brief        peer routine of
  !!                            sp_enefunc_gromacs_mod::define_enefunc_gromacs()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs_peer(grotop,       &
                                         ene_info,     &
                                         domain_index, &
                                         constraints,  &
                                         domain,       &
                                         enefunc)

    ! formal arguments
    type(s_grotop),       intent(in)    :: grotop
    type(s_ene_info),     intent(in)    :: ene_info
    type(s_domain_index), intent(in)    :: domain_index
    type(s_constraints),  intent(inout) :: constraints
    type(s_domain),       intent(inout) :: domain
    type(s_enefunc),      intent(inout) :: enefunc

    ! local variables
    integer               :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncRBDihe,   ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)


    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond_G_peer(grotop, domain_index, domain, enefunc)


      ! angle
      !
      call setup_enefunc_angl_G_peer(grotop, domain_index, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint_G_peer( &
                           grotop, constraints, domain_index, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint_G_peer( &
                           grotop, constraints, domain_index, domain, enefunc)

    end if


    ! dihedral
    !
    call setup_enefunc_dihe_G_peer(grotop, domain_index, domain, enefunc)


    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe_G_peer(grotop, domain_index, domain, enefunc)


    ! improper
    !
    call setup_enefunc_impr_G_peer(grotop, domain_index, domain, enefunc)


    ! nonbonded
    !
    call setup_enefunc_nonb_G_peer(ene_info, grotop, constraints,  &
                                   domain, enefunc)


    ! restraint
    !
    call setup_enefunc_restraints_peer(domain_index, domain, enefunc)


    ! restraints (gromacs)
    !
    ! call setup_enefunc_gro_restraints() ...skip

    return

  end subroutine define_enefunc_gromacs_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_bond()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_G_peer(grotop, domain_index, domain, enefunc)

    ! paramters
    integer, parameter :: WaterIdx(2,3) = reshape((/1,2,1,3,2,3/),shape=(/2,3/))

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3)
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx1, idx2, icel1, icel2, icel_local
    integer                  :: nbond_c

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond_G_Peer> Too many bonds.') 

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                dist (  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%b0 * 10.0_wp
              end if
  
            end if
  
          end do

        else

          water_dist = (/gromol%settles%doh * 10.0_wp, &
                         gromol%settles%doh * 10.0_wp, &
                         gromol%settles%dhh * 10.0_wp/)

          do k = 1, 3

            idx1 = WaterIdx(1,k) + ioffset
            idx2 = WaterIdx(2,k) + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond_G_Peer> Too many bonds.') 

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = 0.0_wp
                dist (  bond(icel_local),icel_local) = water_dist(k)
              end if

            end if

          end do

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_bond_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_angl()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_G_peer(grotop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx1, idx2, idx3, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (angle(icel_local)+1 > MaxAngle) &
                  call error_msg('Setup_Enefunc_Angl_G_Peer> Too many angles.') 

                angle(icel_local) = angle(icel_local) + 1
                list (1:3,angle(icel_local),icel_local) = (/idx1, idx2, idx3/)
                force(    angle(icel_local),icel_local) = &
                                    gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                theta(    angle(icel_local),icel_local) = &
                                    gromol%angls(k)%theta_0 * RAD
              end if

            end if

          end do

        else

          idx1 = 2 + ioffset
          idx2 = 1 + ioffset
          idx3 = 3 + ioffset

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx3)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              if (angle(icel_local)+1 > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl_G_Peer> Too many angles.') 

              angle(icel_local) = angle(icel_local) + 1
              list (1:3,angle(icel_local),icel_local) = (/idx1, idx2, idx3/)
              force(    angle(icel_local),icel_local) = 0.0_wp
              theta(    angle(icel_local),icel_local) = 0.0_wp
            end if

          end if

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_angl_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_dihe()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_G_peer(grotop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types
    notation  = 100

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
              .and. gromol%dihes(k)%func /= 9) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              ndihe => dihe(icel_local)
              ndihe = ndihe + 1

              if (ndihe > MaxDihe) &
              call error_msg('Setup_Enefunc_Dihe_G_Peer> Too many dihedrals.') 

              list (1:4,ndihe,icel_local) = (/idx1, idx2, idx3, idx4/)
              force (   ndihe,icel_local) = gromol%dihes(k)%kp * JOU2CAL
              phase (   ndihe,icel_local) = gromol%dihes(k)%ps * RAD
              period(   ndihe,icel_local) = gromol%dihes(k)%multiplicity
              if (period(ndihe,icel_local) >  enefunc%notation_14types) &
              call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_dihe_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe_G_peer
  !> @brief        peer routine of 
  !!                           sp_enefunc_gromacs_mod::setup_enefunc_rb_dihe()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe_G_peer(grotop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: c(:,:,:)
    integer,            pointer :: dihe(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_rb_dihedral
    list      => enefunc%rb_dihe_list
    c         => enefunc%rb_dihe_c

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 3) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              ndihe => dihe(icel_local)
              ndihe = ndihe + 1

              if (ndihe > MaxDihe) &
            call error_msg('Setup_Enefunc_RB_Dihe_G_Peer> Too many dihedrals.') 

              list (1:4,ndihe,icel_local) = (/idx1, idx2, idx3, idx4/)
              c    (1:6,ndihe,icel_local) = gromol%dihes(k)%c(1:6) * JOU2CAL
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_rb_dihe_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_impr()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_G_peer(grotop, domain_index, domain, enefunc)

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 2) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nimpr => impr(icel_local)
              nimpr = nimpr + 1

              if (nimpr > MaxImpr) &
              call error_msg('Setup_Enefunc_Impr_G_Peer> Too many impropers.') 

              list (1:4,nimpr,icel_local) = (/idx1, idx2, idx3, idx4/)
              force (   nimpr,icel_local) = gromol%dihes(k)%kp * JOU2CAL
              phase (   nimpr,icel_local) = gromol%dihes(k)%ps * RAD
              period(   nimpr,icel_local) = gromol%dihes(k)%multiplicity
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_impr_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::
  !!                                          setup_enefunc_bond_constraint()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint_G_peer(grotop, constraints, &
                                                 domain_index, domain,  enefunc)

    ! paramters
    integer, parameter :: WaterIdx(2,3) = reshape((/1,2,1,3,2,3/),shape=(/2,3/))

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_constraints),  target, intent(inout) :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3)
    integer                  :: i, j, k, m, n
    integer                  :: ioffset
    integer                  :: icel, connect, ih, ih1, ih2, nbond_c
    integer                  :: i1, i2, idx1, idx2, icel1, icel2, icel_local
    character(4)             :: atm1, atm2
    character(6)             :: res1, res2

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:)
    real(wp),           pointer :: HGr_bond_dist(:,:,:,:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:), id_l2g(:,:)
    integer,            pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)


    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist
    connect       =  constraints%connect

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min

    ioffset   = 0
    nbond_c   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds

            i1 = gromol%bonds(k)%atom_idx1
            i2 = gromol%bonds(k)%atom_idx2

            atm1 = gromol%atoms(i1)%atom_type
            atm2 = gromol%atoms(i2)%atom_type

            idx1 = i1 + ioffset
            idx2 = i2 + ioffset

            if (atm1(1:1) /= 'H' .and. atm2(1:1) /= 'H' .and. &
                atm1(1:1) /= 'h' .and. atm2(1:1) /= 'h') then

              icel1 = id_g2l(1,idx1)
              icel2 = id_g2l(1,idx2)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel_local = cell_pair(icel1,icel2)

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (bond(icel_local)+1 > MaxBond) &
        call error_msg('Setup_Enefunc_Bond_Constraint_G_Peer> Too many bonds.') 

                  bond(icel_local) = bond(icel_local) + 1
                  list (1,bond(icel_local),icel_local) = idx1
                  list (2,bond(icel_local),icel_local) = idx2
                  force(  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                  dist (  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%b0 * 10.0_wp
                end if

              end if

            else

              do icel = 1, ncel
                do m = 1, connect

                  do n = 1, HGr_local(m,icel)
                    ih1 = id_l2g(HGr_bond_list(1,n,m,icel),icel)
                    do ih = 1, m
                      ih2 = id_l2g(HGr_bond_list(ih+1,n,m,icel),icel)

                      if (ih1 == idx1 .and. ih2 == idx2 .or. &
                          ih2 == idx1 .and. ih1 == idx2) then

                        nbond_c = nbond_c + 1
                        HGr_bond_dist(ih+1,n,m,icel) = &
                                   gromol%bonds(k)%b0 * 10.0_wp

                        goto 1

                      end if

                    end do
                  end do

                end do
              end do
  1           continue

            end if

          end do

        else

          water_dist = (/gromol%settles%doh * 10.0_wp, &
                         gromol%settles%doh * 10.0_wp, &
                         gromol%settles%dhh * 10.0_wp/)

          do k = 1, 3

            idx1 = WaterIdx(1,k) + ioffset
            idx2 = WaterIdx(2,k) + ioffset

            do icel = 1, ncel
              do m = 1, connect

                do n = 1, HGr_local(m,icel)
                  ih1 = id_l2g(HGr_bond_list(1,n,m,icel),icel)
                  do ih = 1, m
                    ih2 = id_l2g(HGr_bond_list(ih+1,n,m,icel),icel)

                    if (ih1 == idx1 .and. ih2 == idx2 .or. &
                        ih2 == idx1 .and. ih1 == idx2) then

                      HGr_bond_dist(ih+1,n,m,icel) = water_dist(k)
  
                      goto 2

                    end if

                  end do
                end do

              end do
            end do
2           continue

          end do

        end if

        ioffset   = ioffset   + gromol%num_atoms

      end do
    end do
    constraints%num_bonds = nbond_c

    return

  end subroutine setup_enefunc_bond_constraint_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::
  !!                                          setup_enefunc_angl_constraint()
  !! @authors      NT, JJ
  !! @param[in]    
  !! @date         2012/10/24 (JJ) - removing angle list of water molecules
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint_G_peer(grotop, constraints, &
                                                  domain_index, domain, enefunc)

    ! formal arguments
    type(s_grotop),               intent(in)    :: grotop
    type(s_constraints),          intent(in)    :: constraints
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, icel_local
    integer                  :: i1, i2, i3, idx1, idx2, idx3, icel1, icel2
    character(6)             :: res1, res2, res3

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            i1 = gromol%angls(k)%atom_idx1
            i2 = gromol%angls(k)%atom_idx2
            i3 = gromol%angls(k)%atom_idx3

            res1 = gromol%atoms(i1)%residue_name
            res2 = gromol%atoms(i2)%residue_name
            res3 = gromol%atoms(i3)%residue_name

            if (res1(1:4) /= constraints%water_model .and. &
                res2(1:4) /= constraints%water_model .and. &
                res3(1:4) /= constraints%water_model) then

              idx1 = i1 + ioffset
              idx2 = i2 + ioffset
              idx3 = i3 + ioffset

              icel1 = id_g2l(1,idx1)
              icel2 = id_g2l(1,idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel_local = cell_pair(icel1,icel2)

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (angle(icel_local)+1 > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl_Constraint_G_Peer> Too many angles.')

                  angle(icel_local) = angle(icel_local) + 1
                  list (1:3,angle(icel_local),icel_local) = (/idx1,idx2,idx3/)
                  force(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                  theta(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%theta_0 * RAD
                end if

              end if

            end if

          end do

        end if

        ioffset   = ioffset   + gromol%num_atoms

      end do
    end do

    return

  end subroutine setup_enefunc_angl_constraint_G_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_restraints_peer
  !> @brief        peer routine of 
  !                        sp_enefunc_restraints_mod::setup_enefunc_restraints()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_restraints_peer(domain_index, domain, enefunc)

    ! formal arguments
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc


    if (.not. enefunc%restraint) return


    call setup_enefunc_rest_domain_peer(domain_index, domain, enefunc)

    return

  end subroutine setup_enefunc_restraints_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rest_domain_peer
  !> @brief        peer routine of 
  !                       sp_enefunc_restraints_mod::setup_enefunc_rest_domain()
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rest_domain_peer(domain_index, domain, enefunc)

    ! formal arguments
    type(s_domain_index), target, intent(in)    :: domain_index
    type(s_domain),       target, intent(in)    :: domain
    type(s_enefunc),      target, intent(inout) :: enefunc

    ! local variable
    integer                       :: i, k, ig, iatm, icel
    integer                       :: iselatoms, nselatoms_lst

    real(wp),             pointer :: const(:,:)
    real(wp),             pointer :: restraint_force(:,:,:)
    real(wp),             pointer :: restraint_coord(:,:,:)
    integer,              pointer :: id_g2l(:,:)
    integer,              pointer :: ncel, num_funcs
    integer,              pointer :: grouplist(:,:), bondslist(:,:)
    integer,              pointer :: nrestraint(:), restraint_atom(:,:)
    integer,              pointer :: iselatoms_lst(:)


    ncel            => domain%num_cell_local
    id_g2l          => domain%id_g2l

    num_funcs       => enefunc%num_restraintfuncs
    grouplist       => enefunc%restraint_grouplist
    bondslist       => enefunc%restraint_bondslist
    const           => enefunc%restraint_const
    nrestraint      => enefunc%num_restraint

    ! check the size for distance restraint
    !
    k = 0
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncDIST     .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDISTCOM  .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLE    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncANGLECOM .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHED    .or. &
          enefunc%restraint_kind(i) == RestraintsFuncDIHEDCOM ) then

        do ig = 1, enefunc%restraint_funcgrp(i)

          iselatoms_lst => domain_index%selatoms_list(grouplist(ig,i))%i
          nselatoms_lst =  domain_index%selatoms_list(grouplist(ig,i))%n

          k = k + nselatoms_lst

        end do

      end if

    end do
    enefunc%num_atoms_bonds_restraint = k


    call alloc_enefunc(enefunc, EneFuncRest, ncel, 1)

    restraint_atom  => enefunc%restraint_atom
    restraint_force => enefunc%restraint_force
    restraint_coord => enefunc%restraint_coord

    ! position restraint
    !
    do i = 1, num_funcs

      if (enefunc%restraint_kind(i) == RestraintsFuncPOSI) then

        iselatoms_lst => domain_index%selatoms_list(grouplist(1,i))%i
        nselatoms_lst =  domain_index%selatoms_list(grouplist(1,i))%n

        do iselatoms = 1, nselatoms_lst

          iatm = hm_selatoms(iselatoms_lst(iselatoms), grouplist(1,i))
          icel = id_g2l(1, iatm)
          if (icel > 0 .and. icel <= ncel) then
            nrestraint(icel) = nrestraint(icel) + 1
            restraint_atom (    nrestraint(icel),icel) = iatm
            restraint_force(1:4,nrestraint(icel),icel) = const(1:4,i)
            restraint_coord(1,nrestraint(icel),icel) = hm_atom_refcoord(1,iatm)
            restraint_coord(2,nrestraint(icel),icel) = hm_atom_refcoord(2,iatm)
            restraint_coord(3,nrestraint(icel),icel) = hm_atom_refcoord(3,iatm)
          end if

        end do

      end if

    end do


    ! distance restraint
    !

    ! ..skip

    return

  end subroutine setup_enefunc_rest_domain_peer


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_G_peer
  !> @brief        peer routine of sp_enefunc_gromacs_mod::setup_enefunc_nonb()
  !! @authors      NT (modified by JJ)
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_G_peer(ene_info, grotop, constraints, &
                                       domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig, ei, ej, si, sj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, ncel, i, j, k, excl_level
    integer                  :: ix, jx, kx
    integer                  :: cls_local

    integer,    allocatable  :: check_cls(:)
    integer,    allocatable  :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp),   allocatable  :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp),   allocatable  :: nonb_lj6(:,:), nonb_lj12(:,:)

    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ncel                 = domain%num_cell_local
    ELECOEF              = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    allocate(check_cls(nnonb),             &
             atmcls_map_g2l(nnonb),        &
             atmcls_map_l2g(nnonb),        &
             nb14_lj6      (nnonb, nnonb), &
             nb14_lj12     (nnonb, nnonb), &
             nonb_lj6      (nnonb, nnonb), &
             nonb_lj12     (nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb

      do j = 1, nnonb

        ! combination rule
        !
        if (grotop%num_nbonparms == 0) then

          vi = grotop%atomtypes(i)%v
          vj = grotop%atomtypes(j)%v
          wi = grotop%atomtypes(i)%w
          wj = grotop%atomtypes(j)%w

          if (grotop%defaults%combi_rule == 2) then

            si = vi * 10.0_wp
            sj = vj * 10.0_wp

            ei = wi * JOU2CAL
            ej = wj * JOU2CAL

            sig = (si + sj) * 0.5_wp
            eps = sqrt(ei * ej)

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule == 1 or 3

            c6i  = vi * 1000000.0_wp * JOU2CAL
            c6j  = vj * 1000000.0_wp * JOU2CAL

            c12i = wi * 1000000.0_wp * 1000000.0_wp * JOU2CAL
            c12j = wj * 1000000.0_wp * 1000000.0_wp * JOU2CAL

            c6  = sqrt(c6i  * c6j)
            c12 = sqrt(c12i * c12j)

          end if

        else

          vij = 0.0_wp
          wij = 0.0_wp

          do k = 1, grotop%num_nbonparms
            if (grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type2 .or.  &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type2) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (vij == 0.0_wp .or. wij == 0.0_wp) &
       call error_msg('Setup_Enefunc_Nonb_G_Peer_0> combination is not found.')

          if (grotop%defaults%combi_rule == 2) then

            sig = vij * 10.0_wp
            eps = wij * JOU2CAL

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule = 1 or 3

            c6  = vij * 1000000.0_wp * JOU2CAL
            c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL

          end if

        end if

        ! set parameters
        !
        nb14_lj12(i,j) = c12
        nb14_lj6 (i,j) = c6

        nonb_lj12(i,j) = c12
        nonb_lj6 (i,j) = c6

      end do
    end do

    ! check # of exclusion level
    !
    if (enefunc%forcefield == ForcefieldGROMARTINI) then

      !TODO

      enefunc%excl_level = -1

      do i = 1, grotop%num_molss
        excl_level = grotop%molss(i)%moltype%exclude_nbon
        if (enefunc%excl_level == -1) then
          enefunc%excl_level = excl_level
        else if (enefunc%excl_level /= excl_level) then
          call error_msg( &
     'Setup_Enefunc_Nonb_G_Peer_0> multiple "exclude_nbon" is not supported.')
        end if
      end do
    end if

    ! check the usage of atom class
    !
    do i = 1, hm_num_atoms
      k = hm_atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(hm_atom_cls_name(i))//'"')
      endif
      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nnonb
      if (check_cls(i) == 1) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do
    if (constraints%rigid_bond .or. enefunc%table%water_table) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(1:3) = &
        atmcls_map_g2l(domain%water%atom_cls_no(1:3))
      enefunc%table%atom_cls_no_O = atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H = atmcls_map_g2l(enefunc%table%atom_cls_no_H)
    endif

    deallocate(check_cls,     &
               atmcls_map_g2l,&
               atmcls_map_l2g,&
               nb14_lj6,      &
               nb14_lj12,     &
               nonb_lj6,      &
               nonb_lj12)

    enefunc%num_atom_cls = cls_local

    return

  end subroutine setup_enefunc_nonb_G_peer


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_id_g2l
  !> @brief        allocate domain%id_g2l(:,:)
  !! @authors      NT
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_id_g2l(domain, natom)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: natom


    if (allocated(domain%id_g2l)) then
      if (size(domain%id_g2l(1,:)) /= natom) &
        deallocate(domain%id_g2l)
    end if

    if (.not. allocated(domain%id_g2l)) &
      allocate(domain%id_g2l(2, natom))

    domain%id_g2l(1:2,1:natom) = 0

    return

  end subroutine alloc_id_g2l

end module pr_setup_spdyn_peer_mod
