!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_constraints_mod
!> @brief   constraints module
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_constraints_mod

  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_fep_utils_mod

  implicit none
  private

  ! structures
  type, public :: s_cons_info
    logical                  :: rigid_bond      = .false.
    logical                  :: fast_bond       = .false.
    logical                  :: fast_water      = .false.
    integer                  :: hydrogen_type   = ConstraintAtomName
    integer                  :: shake_iteration = 500
    real(wp)                 :: shake_tolerance = 1.0e-10_dp
    integer                  :: lincs_iteration = 1
    integer                  :: lincs_order     = 4
    character(5)             :: water_model     = 'TIP3'
    real(wp)                 :: hydrogen_mass_upper_bound = 2.1_wp
  end type s_cons_info

  ! parameters
  integer, public, parameter :: ConstraintModeLEAP  = 1
  integer, public, parameter :: ConstraintModeVVER1 = 2
  integer, public, parameter :: ConstraintModeVVER2 = 3

  ! subroutines
  public  :: show_ctrl_constraints
  public  :: read_ctrl_constraints
  public  :: setup_constraints
  public  :: setup_constraints_pio
  public  :: compute_constraints
  public  :: water_force_redistribution
  public  :: decide_dummy
  public  :: compute_settle_min
  private :: setup_fast_water
  private :: setup_fast_water_pio
  private :: setup_rigid_bond
  private :: setup_rigid_bond_pio
  private :: compute_settle
  private :: compute_shake
  private :: compute_rattle_fast_vv1
  private :: compute_rattle_vv1
  private :: compute_rattle_fast_vv2
  private :: compute_rattle_vv2

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_constraints
  !> @brief        show CONSTRAINTS section usage
  !! @authors      JJ
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_constraints(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO         # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') '# shake_iteration = 500      # max number of SHAKE/RATTLE iterations'
        write(MsgOut,'(A)') '# shake_tolerance = 1.0e-8   # SHAKE/RATTLE tolerance (Ang)'
        write(MsgOut,'(A)') '# water_model     = TIP3    # water model'
        write(MsgOut,'(A)') '# hydrogen_mass_upper_bound  = 2.1    # upper mass limit to define the hydrogen atom'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[CONSTRAINTS]'
        write(MsgOut,'(A)') 'rigid_bond    = NO         # constraints all bonds involving hydrogen'
        write(MsgOut,'(A)') ' '

      end select

    end if


    return

  end subroutine show_ctrl_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_constraint
  !> @brief        read CONSTRAINTS section in the control file
  !! @authors      JJ
  !! @param[in]    handle    :unit number
  !! @param[out]   cons_info :CONSTRAINTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_constraints(handle, cons_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Constraints'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_cons_info),       intent(inout) :: cons_info
  

    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'rigid_bond',      &
                               cons_info%rigid_bond)
    call read_ctrlfile_logical(handle, Section, 'fast_water',      &
                               cons_info%fast_water)
    call read_ctrlfile_integer(handle, Section, 'shake_iteration', &
                               cons_info%shake_iteration)
    call read_ctrlfile_real   (handle, Section, 'shake_tolerance', &
                               cons_info%shake_tolerance)
    call read_ctrlfile_string (handle, Section, 'water_model',     &
                               cons_info%water_model)
    call read_ctrlfile_type(handle, Section, 'hydrogen_type',      &
                               cons_info%hydrogen_type, ConstraintAtomType)
    call read_ctrlfile_real   (handle, Section, 'hydrogen_mass_upper_bound', &
                               cons_info%hydrogen_mass_upper_bound)

    call end_ctrlfile_section(handle)


    ! fast_water is true if rigid_bond is true
    !
    if (cons_info%rigid_bond) cons_info%fast_water = .true.

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Constraints> Parameters for Constraints'

      if (cons_info%rigid_bond) then
        write(MsgOut,'(A30)')                                    &
              '  rigid_bond      =        yes'
        write(MsgOut,'(A20,I10,A20,E10.3)')                      &
              '  shake_iteration = ', cons_info%shake_iteration, &
              '  shake_tolerance = ', cons_info%shake_tolerance

        if (cons_info%fast_bond) then
          write(MsgOut,'(A30)')                                  &
                '  fast_bond       =        yes'
          write(MsgOut,'(A20,I10,A20,I10)')                        &
                '  lincs_iteration = ', cons_info%lincs_iteration, &
                '  lincs_order     = ', cons_info%lincs_order

        else
          write(MsgOut,'(A30)')                                  &
                '  fast_bond       =         no'
        end if

        if (cons_info%fast_water) then
          write(MsgOut,'(A30,A20,A10)')                          &
                '  fast_water      =        yes',                &
                '  water_model     = ', trim(cons_info%water_model)
        else
          write(MsgOut,'(A30)')                                  &
                '  fast_water      =         no'
        end if

        if (cons_info%hydrogen_type==ConstraintAtomName) then
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =       name'
        else if (cons_info%hydrogen_type==ConstraintAtomMass) then
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =       mass'
        else 
          write(MsgOut,'(A30)')                                    &
                '  hydrogen_type   =  name|mass'
        end if

      else
        write(MsgOut,'(A30)') '  rigid_bond      =         no'
      end if

      write(MsgOut,'(A)') ' '

    end if


    ! error check
    !
    if (main_rank) then

      if (cons_info%fast_bond .and. &
          .not. cons_info%rigid_bond) then
        call error_msg( &
          'Read_Ctrl_Constraints> rigid_bond must be YES, if fast_bond is YES')
      end if

!TODO
!      if (cons_info%fast_water .and. &
!          .not. cons_info%rigid_bond) then
!        call error_msg( &
!          'Read_Ctrl_Constraints> rigid_bond must be YES, if fast_water is YES')
!      end if

    end if


    return    

  end subroutine read_ctrl_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints
  !> @brief        setup constraints for domain decomposition
  !! @authors      JJ
  !! @param[in]    cons_info :CONSTRAINTS section control parameters information
  !! @param[in]    par         : PAR information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints(cons_info, par, prmtop, grotop, &
                               molecule, enefunc, constraints)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints


    constraints%rigid_bond      = cons_info%rigid_bond
    constraints%fast_bond       = cons_info%fast_bond
    constraints%fast_water      = cons_info%fast_water
    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%lincs_order     = cons_info%lincs_order
    constraints%water_model     = cons_info%water_model
    constraints%hydrogen_type   = cons_info%hydrogen_type

    if (constraints%rigid_bond) then

      if (constraints%hydrogen_type == ConstraintAtomName  &
         .and. molecule%special_hydrogen) &
      call error_msg('Setup_Constraints> Non ordinary hydrogen name is not'//&
                         ' allowed. If you want, use hydrogen_type option.')

      ! setup SETTLE
      !
      if (constraints%fast_water .and. enefunc%table%num_water > 0) then

        if (constraints%tip4) then

          call setup_fast_water_tip4(par, prmtop, grotop, &
                                     molecule, enefunc, constraints)
        else

          call setup_fast_water(par, prmtop, grotop, &
                                molecule, enefunc, constraints)

        end if

        ! update number of degrees of freedom
        !
        if (constraints%tip4) then

          call update_num_deg_freedom('After setup of SETTLE',    &
                                      -6*enefunc%table%num_water, &
                                      molecule%num_deg_freedom)

        else

          call update_num_deg_freedom('After setup of SETTLE',    &
                                      -3*enefunc%table%num_water, &
                                      molecule%num_deg_freedom)

        end if

      end if

      ! setup SHAKE and RATTLE
      !
      call setup_rigid_bond(par, molecule, enefunc, constraints)

      ! update number of degrees of freedom
      !
      if (constraints%num_bonds > 0) then
        call update_num_deg_freedom('After setup of SHAKE/RATTLE', &
                                    -constraints%num_bonds,        &
                                    molecule%num_deg_freedom)

      end if

    else if (constraints%fast_water .and. .not.enefunc%table%tip4 &
             .and. enefunc%table%num_water > 0) then

      call setup_fast_water(par, prmtop, grotop, &
                            molecule, enefunc, constraints)

    else if (enefunc%table%tip4 .and. &
             enefunc%table%num_water > 0) then

      call setup_fast_water_tip4(par, prmtop, grotop, &
                                 molecule, enefunc, constraints)

    end if
   
    return

  end subroutine setup_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_constraints_pio
  !> @brief        setup constraints for domain decomposition
  !! @authors      JJ
  !! @param[in]    cons_info :CONSTRAINTS section control parameters information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_constraints_pio(cons_info, pio_restart, enefunc,  &
                                   constraints, domain)

    ! formal arguments
    type(s_cons_info),       intent(in)    :: cons_info
    logical,                 intent(in)    :: pio_restart
    type(s_constraints),     intent(inout) :: constraints
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_domain),          intent(inout) :: domain


    constraints%rigid_bond      = cons_info%rigid_bond
    constraints%fast_bond       = cons_info%fast_bond
    constraints%fast_water      = cons_info%fast_water
    constraints%shake_iteration = cons_info%shake_iteration
    constraints%shake_tolerance = cons_info%shake_tolerance
    constraints%lincs_iteration = cons_info%lincs_iteration
    constraints%lincs_order     = cons_info%lincs_order
    constraints%water_model     = cons_info%water_model

    if (constraints%rigid_bond) then

      ! setup SETTLE
      !
      if (constraints%fast_water) then

        call setup_fast_water_pio(constraints)

        ! update number of degrees of freedom
        !
        call update_num_deg_freedom('After setup of SETTLE',    &
                                    -3*enefunc%table%num_water, &
                                    domain%num_deg_freedom)

      end if

      ! setup SHAKE and RATTLE
      !
      call setup_rigid_bond_pio(constraints)

      ! update number of degrees of freedom
      !
      if (constraints%num_bonds > 0) then
        call update_num_deg_freedom('After setup of SHAKE/RATTLE', &
                                    -constraints%num_bonds,        &
                                    domain%num_deg_freedom)

      end if

    end if
   
    return

  end subroutine setup_constraints_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_constraints
  !> @brief        update coordinates according to constraints
  !! @authors      JJ
  !! @param[in]    cons_mode   : constraint mode
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[inout] coord_ref   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_constraints(cons_mode, vel_update, dt, coord_ref, &
                                 domain, constraints, coord, vel, virial)

    ! formal arguments
    integer,                 intent(in)    :: cons_mode
    logical,                 intent(in)    :: vel_update
    real(dp),                intent(in)    :: dt
    real(dp),                intent(inout) :: coord_ref(:,:,:)
    type(s_domain),          intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(dp),                intent(inout) :: coord(:,:,:)
    real(dp),                intent(inout) :: vel(:,:,:)
    real(dp),                intent(inout) :: virial(:,:)


    call timer(TimerConstraint, TimerOn)

    ! constraints
    !
    select case (cons_mode)

    case (ConstraintModeLEAP)

      virial(1:3,1:3) = 0.0_dp

      call compute_settle(vel_update, dt, coord_ref, domain, &
                          constraints, coord, vel, virial)

      call compute_shake (vel_update, dt, coord_ref, domain, &
                          constraints, coord, vel, virial)

      if (constraints%tip4) &
        call decide_dummy(domain, constraints, coord)

      virial(1:3,1:3) = virial(1:3,1:3)/dt**2

    case (ConstraintModeVVER1)

      call compute_rattle_fast_vv1(dt, coord_ref, &
                                   domain, constraints, coord, vel)

      call compute_rattle_vv1     (dt, coord_ref, &
                                   domain, constraints, coord, vel)

      if (constraints%tip4) &
        call decide_dummy(domain, constraints, coord)

    case (ConstraintModeVVER2)

      call compute_rattle_fast_vv2(domain, constraints, coord, vel)

      call compute_rattle_vv2     (domain, constraints, coord, vel)

    end select

    ! FEP: synchronize single B with single A
    if (domain%fep_use) then
      call sync_single_fep(domain, coord)
      call sync_single_fep(domain, vel)
    end if

    call timer(TimerConstraint, TimerOff)

    return

  end subroutine compute_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water
  !> @brief        setup parameter for bond in water
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water(par, prmtop, grotop, &
                              molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
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
    list(1:3) = enefunc%table%water_list(1:3,1)

    ! charmm
    if (par%num_bonds > 0) then

      ci1 = molecule%atom_cls_name(list(1))
      ci2 = molecule%atom_cls_name(list(2))
      ci3 = molecule%atom_cls_name(list(3))

      do i = 1, par%num_bonds
        if (ci1 == par%bond_atom_cls(1,i) .and. &
            ci2 == par%bond_atom_cls(2,i) .or.  &
            ci1 == par%bond_atom_cls(2,i) .and. &
            ci2 == par%bond_atom_cls(1,i) ) then

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

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_dp
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_dp
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

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_dp
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rHH = gromol%settles%dhh * 10.0_dp
          goto 2

        end if

      end do
2     continue

    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Fast_Water> Setup constraints for SETTLE'
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(O-H)         = ', constraints%water_rOH,  &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(H-H)         = ', constraints%water_rHH,  &
           '  mass(H)         = ', constraints%water_massH
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_tip4
  !> @brief        setup parameter for bond in water
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_tip4(par, prmtop, grotop, &
                                   molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, i1, i2, list(4)
    integer                  :: ioffset, j, k
    real(wp)                 :: a, b
    character(6)             :: ci1, ci2, ci3, ci4

    type(s_grotop_mol), pointer :: gromol

    ! mass
    !
    constraints%water_massO = enefunc%table%mass_O
    constraints%water_massH = enefunc%table%mass_H

    ! min distance
    !
    list(1:4) = enefunc%table%water_list(1:4,1)

    ! charmm
    !
    if (par%num_bonds > 0) then

      ci1 = molecule%atom_cls_name(list(1))
      ci2 = molecule%atom_cls_name(list(2))
      ci3 = molecule%atom_cls_name(list(3))
      ci4 = molecule%atom_cls_name(list(4))

      do i = 1, par%num_bonds
        if (ci1 == par%bond_atom_cls(1,i) .and. &
            ci2 == par%bond_atom_cls(2,i) .or.  &
            ci1 == par%bond_atom_cls(2,i) .and. &
            ci2 == par%bond_atom_cls(1,i) ) then

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
        if (ci1 == par%bond_atom_cls(1,i) .and. &
            ci4 == par%bond_atom_cls(2,i) .or.  &
            ci1 == par%bond_atom_cls(2,i) .and. &
            ci4 == par%bond_atom_cls(1,i) ) then

          constraints%water_rOD = par%bond_dist_min(i)
          exit

        end if
      end do

    ! amber
    !
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

        if (list(1) == i1 .and. list(4) == i2 .or. &
            list(4) == I1 .and. list(1) == i2) then

          constraints%water_rOD = &
              prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))
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

                constraints%water_rOH = gromol%bonds(k)%b0 * 10.0_dp
                goto 1

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else

          constraints%water_rOH = gromol%settles%doh * 10.0_dp
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

                constraints%water_rHH = gromol%bonds(k)%b0 * 10.0_dp
                goto 2

              end if
            end do

            ioffset = ioffset + gromol%num_atoms

          end do

        else


          constraints%water_rHH = gromol%settles%dhh * 10.0_dp
          goto 2

        end if

      end do

2     ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol

        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_vsites3
            i1 = gromol%vsites3(k)%func

            if (i1 == 1) then

              ioffset = 1
              a = gromol%vsites3(k)%a
              b = gromol%vsites3(k)%b
              if (a /= b) call error_msg('TIP4P molecule is not symmetric')
              constraints%water_rOD = 2.0_wp * a  &
                                * sqrt(constraints%water_rOH**2 &
                                       -0.25_wp*constraints%water_rHH**2)
            end if
            if (ioffset == 1) exit
          end do

          if (ioffset == 1) exit

        end do

        if (ioffset == 1) exit
      end do

      if (ioffset /= 1) call error_msg('Virtual site should be defined &
                                     & when using TIP4P')

    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Fast_Water> Setup constraints for SETTLE'
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(O-H)         = ', constraints%water_rOH,  &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4)')                          &
           '  r0(O-D)         = ', constraints%water_rOD
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(H-H)         = ', constraints%water_rHH,  &
           '  mass(H)         = ', constraints%water_massH
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water_tip4

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fast_water_pio
  !> @brief        setup parameter for bond in water
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fast_water_pio(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Fast_Water_Pio> Constraints for SETTLE'
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(O-H)         = ', constraints%water_rOH,  &
           '  mass(O)         = ', constraints%water_massO
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                &
           '  r0(H-H)         = ', constraints%water_rHH,  &
           '  mass(H)         = ', constraints%water_massH
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_fast_water_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigid_bond
  !> @brief        setup rigid bonds for SHAKE and RATTLE in domain
  !!               decomposition
  !! @authors      JJ
  !! @param[in]    par         : PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigid_bond(par, molecule, enefunc, constraints)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_constraints),     intent(inout) :: constraints


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Rigid_Bond> Setup constrains for SHAKE and RATTLE'
      write(MsgOut,'(A20,I10)') &
           '  num_rigid_bonds = ', constraints%num_bonds
      write(MsgOut,'(A)') &
           ' '
    end if

    return

  end subroutine setup_rigid_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigid_bond_pio
  !> @brief        setup rigid bonds for SHAKE and RATTLE in domain
  !!               decomposition
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigid_bond_pio(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Rigid_Bond_Pio> Constrains for SHAKE and RATTLE'
      write(MsgOut,'(A20,I10)') &
           '  num_rigid_bonds = ', constraints%num_bonds
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_rigid_bond_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle(vel_update, dt, coord_old, domain, &
                            constraints, coord, vel, virial)

    ! formal arguments
    logical,                 intent(in)    :: vel_update
    real(dp),                intent(in)    :: dt
    real(dp),                intent(in)    :: coord_old(:,:,:)
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(dp),                intent(inout) :: coord(:,:,:)
    real(dp),                intent(inout) :: vel(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)

    ! local variables
    real(dp)                 :: rHH, rOH, massO, massH
    real(dp)                 :: mass(1:3), mass_H2O, ra, rb, rc, inv_ra, inv_dt
    real(dp)                 :: x0(1:3,1:3), x1(1:3,1:3), x3(1:3,1:3)
    real(dp)                 :: delt(1:3,1:3), rf(1:3,1:3,1:3)
    real(dp)                 :: com0(1:3), com1(1:3)
    real(dp)                 :: oh1(1:3), oh2(1:3)
    real(dp)                 :: Xaxis(1:3), Yaxis(1:3), Zaxis(1:3)
    real(dp)                 :: Xaxis2(1:3), Yaxis2(1:3), Zaxis2(1:3)
    real(dp)                 :: mtrx(1:3,1:3)
    real(dp)                 :: xp0(1:3,1:3), xp1(1:3,1:3), xp2(1:3,1:3)
    real(dp)                 :: xp3(1:3,1:3)
    real(dp)                 :: dxp21(1:2), dxp23(1:2), dxp31(1:2)
    real(dp)                 :: sin_phi, cos_phi, sin_psi, cos_psi, tmp
    real(dp)                 :: rb_cosphi, rb_sinphi
    real(dp)                 :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(dp)                 :: alpha, beta, gamma, al2bt2, sin_theta, cos_theta
    real(dp)                 :: viri_local(1:3,1:3)
    real(dp)                 :: virial_omp(3,3,nthread)
    integer                  :: i, ix, iatom(1:3), i1, j1, k1
    integer                  :: id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mass(1)    =  massO
    mass(2)    =  massH
    mass(3)    =  massH
    mass_H2O   =  massO + 2.0_dp*massH
    rc         =  0.5_dp * rHH 
    ra         =  2.0_dp * dsqrt(rOH*rOH-rc*rc) * massH / mass_H2O
    rb         =  dsqrt(rOH*rOH-rc*rc) - ra
    inv_ra     =  1.0_dp / ra
    inv_dt     =  1.0_dp / dt

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, i1, j1, k1, x0, x1, x3, com0, com1, oh1, oh2,  &
    !$omp         delt, rf, Xaxis, Yaxis, Zaxis, Xaxis2, Yaxis2, Zaxis2,       &
    !$omp         mtrx, xp0, xp1, xp2, xp3, tmp, sin_phi, cos_phi, sin_psi,    &
    !$omp         cos_psi, rb_cosphi, rb_sinphi, rc_sinpsi_sinphi,             &
    !$omp         rc_sinpsi_cosphi, alpha, beta, gamma, dxp21, dxp23, dxp31,   &
    !$omp         al2bt2, sin_theta, cos_theta, viri_local, id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    viri_local(1:3,1:3) = 0.0_dp

    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)
      com0(1:3)   = 0.0_dp
      com1(1:3)  = 0.0_dp

      do k1 = 1, 3
        x0(1:3,k1) = coord_old(1:3,iatom(k1),i)
        x1(1:3,k1) = coord    (1:3,iatom(k1),i)
        com0(1:3)  = com0(1:3) + x0(1:3,k1)*mass(k1)
        com1(1:3)  = com1(1:3) + x1(1:3,k1)*mass(k1)
      end do
      com0(1:3) = com0(1:3) / mass_H2O
      com1(1:3) = com1(1:3) / mass_H2O

      do k1 = 1, 3
        x0(1:3,k1) = x0(1:3,k1) - com0(1:3)
        x1(1:3,k1) = x1(1:3,k1) - com1(1:3)
      end do

      oh1(1:3) = x0(1:3,2) - x0(1:3,1)
      oh2(1:3) = x0(1:3,3) - x0(1:3,1)

      Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
      Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
      Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
      Xaxis(1) = x1(2,1)*Zaxis(3) - x1(3,1)*Zaxis(2)
      Xaxis(2) = x1(3,1)*Zaxis(1) - x1(1,1)*Zaxis(3)
      Xaxis(3) = x1(1,1)*Zaxis(2) - x1(2,1)*Zaxis(1)
      Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
      Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
      Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)

      Xaxis2(1:3) = Xaxis(1:3) * Xaxis(1:3)
      Yaxis2(1:3) = Yaxis(1:3) * Yaxis(1:3)
      Zaxis2(1:3) = Zaxis(1:3) * Zaxis(1:3)
      mtrx(1:3,1) = Xaxis(1:3) / dsqrt(Xaxis2(1)+Xaxis2(2)+Xaxis2(3))
      mtrx(1:3,2) = Yaxis(1:3) / dsqrt(Yaxis2(1)+Yaxis2(2)+Yaxis2(3))
      mtrx(1:3,3) = Zaxis(1:3) / dsqrt(Zaxis2(1)+Zaxis2(2)+Zaxis2(3))

      xp0(1:3,1:3) = 0.0_dp
      xp1(1:3,1:3) = 0.0_dp
      do i1 = 1, 3
        do k1 = 1, 3
          do j1 = 1, 3
            xp0(k1,i1) = xp0(k1,i1) + mtrx(j1,i1)*x0(j1,k1)
            xp1(k1,i1) = xp1(k1,i1) + mtrx(j1,i1)*x1(j1,k1)
          end do
        end do
      end do

      sin_phi = xp1(1,3)*inv_ra
      tmp = 1.0_dp - sin_phi*sin_phi
      if (tmp > 0.0_dp) then
        cos_phi = dsqrt(tmp)
      else
        cos_phi = 0.0_dp
      end if

      sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
      tmp = 1.0_dp - sin_psi*sin_psi
      if (tmp > 0.0_dp) then
        cos_psi = dsqrt(tmp)
      else
        cos_psi = 0.0_dp
      end if

      rb_cosphi = rb * cos_phi
      rb_sinphi = rb * sin_phi
      rc_sinpsi_sinphi = rc * sin_psi * sin_phi
      rc_sinpsi_cosphi = rc * sin_psi * cos_phi

      xp2(1,2) =   ra * cos_phi
      xp2(1,3) =   ra * sin_phi
      xp2(2,1) = - rc * cos_psi
      xp2(2,2) = - rb_cosphi - rc_sinpsi_sinphi
      xp2(2,3) = - rb_sinphi + rc_sinpsi_cosphi
      xp2(3,1) = - xp2(2,1)
      xp2(3,2) = - rb_cosphi + rc_sinpsi_sinphi
      xp2(3,3) = - rb_sinphi - rc_sinpsi_cosphi

      dxp21(1:2) = xp0(2,1:2) - xp0(1,1:2)
      dxp23(1:2) = xp0(2,1:2) - xp0(3,1:2)
      dxp31(1:2) = xp0(3,1:2) - xp0(1,1:2)
      alpha =   xp2(2,1)*dxp23(1) + dxp21(2)*xp2(2,2) + dxp31(2)*xp2(3,2)
      beta  = - xp2(2,1)*dxp23(2) + dxp21(1)*xp2(2,2) + dxp31(1)*xp2(3,2)
      gamma =   dxp21(1) * xp1(2,2) - xp1(2,1) * dxp21(2) &
              + dxp31(1) * xp1(3,2) - xp1(3,1) * dxp31(2)

      al2bt2 = alpha*alpha + beta*beta
      sin_theta = (alpha*gamma - beta*dsqrt(al2bt2 - gamma*gamma))/al2bt2
      tmp = 1.0_dp - sin_theta*sin_theta
      if (tmp > 0.0_dp) then
        cos_theta = dsqrt(tmp)
      else
        cos_theta = 0.0_dp
      end if

      xp3(1,1) = - xp2(1,2)*sin_theta
      xp3(1,2) =   xp2(1,2)*cos_theta
      xp3(1,3) =   xp2(1,3)
      xp3(2,1) =   xp2(2,1)*cos_theta - xp2(2,2)*sin_theta
      xp3(2,2) =   xp2(2,1)*sin_theta + xp2(2,2)*cos_theta
      xp3(2,3) =   xp2(2,3)
      xp3(3,1) =   xp2(3,1)*cos_theta - xp2(3,2)*sin_theta
      xp3(3,2) =   xp2(3,1)*sin_theta + xp2(3,2)*cos_theta
      xp3(3,3) =   xp2(3,3)

      x3(1:3,1:3) = 0.0_dp
      do k1 = 1, 3
        do i1 = 1, 3
          do j1 = 1, 3
            x3(i1,k1) = x3(i1,k1) + mtrx(i1,j1)*xp3(k1,j1)
          end do
        end do
        coord(1:3,iatom(k1),i) = x3(1:3,k1) + com1(1:3)
      end do

      do k1 = 1, 3
        delt(1:3,k1) = x3(1:3,k1) - x1(1:3,k1)
        do i1 = 1, 3
          do j1 = 1, 3
            rf(i1,j1,k1) = - coord_old(i1,iatom(k1),i)*delt(j1,k1)*mass(k1)
          end do
        end do
      end do

      do i1 = 1, 3
        do j1 = 1, 3
          viri_local(1:3,j1) = viri_local(1:3,j1) + rf(1:3,j1,i1)
        end do
      end do

      if (vel_update) then
        do k1 = 1, 3
          vel(1:3,iatom(k1),i) = vel(1:3,iatom(k1),i) + delt(1:3,k1)*inv_dt
        end do
      end if

      end do

    end do

    do k1 = 1, 3
      virial_omp(1:3,k1,id+1) = virial_omp(1:3,k1,id+1) - viri_local(1:3,k1)
    end do

    !$omp end parallel

    do i = 1, nthread
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    return

  end subroutine compute_settle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_settle_min
  !> @brief        SETTLE for three-site water (only coordinaets)
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_settle_min(coord_old, domain, constraints, coord)

    ! formal arguments
    real(dp),                intent(in)    :: coord_old(:,:,:)
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(inout) :: constraints
    real(dp),                intent(inout) :: coord(:,:,:)

    ! local variables
    real(dp)                 :: rHH, rOH, massO, massH
    real(dp)                 :: mass(1:3), mass_H2O, ra, rb, rc, inv_ra, inv_dt
    real(dp)                 :: x0(1:3,1:3), x1(1:3,1:3), x3(1:3,1:3)
    real(dp)                 :: delt(1:3,1:3), rf(1:3,1:3,1:3)
    real(dp)                 :: com0(1:3), com1(1:3)
    real(dp)                 :: oh1(1:3), oh2(1:3)
    real(dp)                 :: Xaxis(1:3), Yaxis(1:3), Zaxis(1:3)
    real(dp)                 :: Xaxis2(1:3), Yaxis2(1:3), Zaxis2(1:3)
    real(dp)                 :: mtrx(1:3,1:3)
    real(dp)                 :: xp0(1:3,1:3), xp1(1:3,1:3), xp2(1:3,1:3)
    real(dp)                 :: xp3(1:3,1:3)
    real(dp)                 :: dxp21(1:2), dxp23(1:2), dxp31(1:2)
    real(dp)                 :: sin_phi, cos_phi, sin_psi, cos_psi, tmp
    real(dp)                 :: rb_cosphi, rb_sinphi
    real(dp)                 :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(dp)                 :: alpha, beta, gamma, al2bt2, sin_theta, cos_theta
    real(dp)                 :: viri_local(1:3,1:3)
    real(dp)                 :: virial_omp(3,3,nthread)
    integer                  :: i, ix, iatom(1:3), i1, j1, k1
    integer                  :: id, omp_get_thread_num

    integer,         pointer :: nwater(:), water_list(:,:,:), ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mass(1)    =  massO
    mass(2)    =  massH
    mass(3)    =  massH
    mass_H2O   =  massO + 2.0_dp*massH
    rc         =  0.5_dp * rHH 
    ra         =  2.0_dp * dsqrt(rOH*rOH-rc*rc) * massH / mass_H2O
    rb         =  dsqrt(rOH*rOH-rc*rc) - ra
    inv_ra     =  1.0_dp / ra

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, i1, j1, k1, x0, x1, x3, com0, com1, oh1, oh2,  &
    !$omp         delt, rf, Xaxis, Yaxis, Zaxis, Xaxis2, Yaxis2, Zaxis2,       &
    !$omp         mtrx, xp0, xp1, xp2, xp3, tmp, sin_phi, cos_phi, sin_psi,    &
    !$omp         cos_psi, rb_cosphi, rb_sinphi, rc_sinpsi_sinphi,             &
    !$omp         rc_sinpsi_cosphi, alpha, beta, gamma, dxp21, dxp23, dxp31,   &
    !$omp         al2bt2, sin_theta, cos_theta, viri_local, id)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)
      com0(1:3)   = 0.0_dp
      com1(1:3)  = 0.0_dp

      do k1 = 1, 3
        x0(1:3,k1) = coord_old(1:3,iatom(k1),i)
        x1(1:3,k1) = coord    (1:3,iatom(k1),i)
        com0(1:3)  = com0(1:3) + x0(1:3,k1)*mass(k1)
        com1(1:3)  = com1(1:3) + x1(1:3,k1)*mass(k1)
      end do
      com0(1:3) = com0(1:3) / mass_H2O
      com1(1:3) = com1(1:3) / mass_H2O

      do k1 = 1, 3
        x0(1:3,k1) = x0(1:3,k1) - com0(1:3)
        x1(1:3,k1) = x1(1:3,k1) - com1(1:3)
      end do

      oh1(1:3) = x0(1:3,2) - x0(1:3,1)
      oh2(1:3) = x0(1:3,3) - x0(1:3,1)

      Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
      Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
      Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
      Xaxis(1) = x1(2,1)*Zaxis(3) - x1(3,1)*Zaxis(2)
      Xaxis(2) = x1(3,1)*Zaxis(1) - x1(1,1)*Zaxis(3)
      Xaxis(3) = x1(1,1)*Zaxis(2) - x1(2,1)*Zaxis(1)
      Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
      Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
      Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)

      Xaxis2(1:3) = Xaxis(1:3) * Xaxis(1:3)
      Yaxis2(1:3) = Yaxis(1:3) * Yaxis(1:3)
      Zaxis2(1:3) = Zaxis(1:3) * Zaxis(1:3)
      mtrx(1:3,1) = Xaxis(1:3) / dsqrt(Xaxis2(1)+Xaxis2(2)+Xaxis2(3))
      mtrx(1:3,2) = Yaxis(1:3) / dsqrt(Yaxis2(1)+Yaxis2(2)+Yaxis2(3))
      mtrx(1:3,3) = Zaxis(1:3) / dsqrt(Zaxis2(1)+Zaxis2(2)+Zaxis2(3))

      xp0(1:3,1:3) = 0.0_dp
      xp1(1:3,1:3) = 0.0_dp
      do i1 = 1, 3
        do k1 = 1, 3
          do j1 = 1, 3
            xp0(k1,i1) = xp0(k1,i1) + mtrx(j1,i1)*x0(j1,k1)
            xp1(k1,i1) = xp1(k1,i1) + mtrx(j1,i1)*x1(j1,k1)
          end do
        end do
      end do

      sin_phi = xp1(1,3)*inv_ra
      tmp = 1.0_dp - sin_phi*sin_phi
      if (tmp > 0.0_dp) then
        cos_phi = dsqrt(tmp)
      else
        cos_phi = 0.0_dp
      end if

      sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
      tmp = 1.0_dp - sin_psi*sin_psi
      if (tmp > 0.0_dp) then
        cos_psi = dsqrt(tmp)
      else
        cos_psi = 0.0_dp
      end if

      rb_cosphi = rb * cos_phi
      rb_sinphi = rb * sin_phi
      rc_sinpsi_sinphi = rc * sin_psi * sin_phi
      rc_sinpsi_cosphi = rc * sin_psi * cos_phi

      xp2(1,2) =   ra * cos_phi
      xp2(1,3) =   ra * sin_phi
      xp2(2,1) = - rc * cos_psi
      xp2(2,2) = - rb_cosphi - rc_sinpsi_sinphi
      xp2(2,3) = - rb_sinphi + rc_sinpsi_cosphi
      xp2(3,1) = - xp2(2,1)
      xp2(3,2) = - rb_cosphi + rc_sinpsi_sinphi
      xp2(3,3) = - rb_sinphi - rc_sinpsi_cosphi

      dxp21(1:2) = xp0(2,1:2) - xp0(1,1:2)
      dxp23(1:2) = xp0(2,1:2) - xp0(3,1:2)
      dxp31(1:2) = xp0(3,1:2) - xp0(1,1:2)
      alpha =   xp2(2,1)*dxp23(1) + dxp21(2)*xp2(2,2) + dxp31(2)*xp2(3,2)
      beta  = - xp2(2,1)*dxp23(2) + dxp21(1)*xp2(2,2) + dxp31(1)*xp2(3,2)
      gamma =   dxp21(1) * xp1(2,2) - xp1(2,1) * dxp21(2) &
              + dxp31(1) * xp1(3,2) - xp1(3,1) * dxp31(2)

      al2bt2 = alpha*alpha + beta*beta
      sin_theta = (alpha*gamma - beta*dsqrt(al2bt2 - gamma*gamma))/al2bt2
      tmp = 1.0_dp - sin_theta*sin_theta
      if (tmp > 0.0_dp) then
        cos_theta = dsqrt(tmp)
      else
        cos_theta = 0.0_dp
      end if

      xp3(1,1) = - xp2(1,2)*sin_theta
      xp3(1,2) =   xp2(1,2)*cos_theta
      xp3(1,3) =   xp2(1,3)
      xp3(2,1) =   xp2(2,1)*cos_theta - xp2(2,2)*sin_theta
      xp3(2,2) =   xp2(2,1)*sin_theta + xp2(2,2)*cos_theta
      xp3(2,3) =   xp2(2,3)
      xp3(3,1) =   xp2(3,1)*cos_theta - xp2(3,2)*sin_theta
      xp3(3,2) =   xp2(3,1)*sin_theta + xp2(3,2)*cos_theta
      xp3(3,3) =   xp2(3,3)

      x3(1:3,1:3) = 0.0_dp
      do k1 = 1, 3
        do i1 = 1, 3
          do j1 = 1, 3
            x3(i1,k1) = x3(i1,k1) + mtrx(i1,j1)*xp3(k1,j1)
          end do
        end do
        coord(1:3,iatom(k1),i) = x3(1:3,k1) + com1(1:3)
      end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_settle_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_shake
  !> @brief        SHAKE for rigid bonds
  !! @authors      JJ
  !! @param[in]    vel_update  : flag for update velocity or not
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !! @param[inout] virial      : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_shake(vel_update, dt, coord_old, domain, &
                            constraints, coord, vel, virial)

    ! formal arguments
    logical,                     intent(in)    :: vel_update
    real(dp),                    intent(in)    :: dt
    real(dp),                    intent(in)    :: coord_old(:,:,:)
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(dp),                    intent(inout) :: coord(:,:,:)
    real(dp),                    intent(inout) :: vel(:,:,:)
    real(dp),                    intent(inout) :: virial(3,3)

    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: x12_old, y12_old, z12_old
    real(dp)                     :: dist2, r
    real(dp)                     :: factor, g12, g12m1, g12m2, v12m1, v12m2
    real(dp)                     :: viri(1:3), fx, fy, fz
    real(dp)                     :: virial_omp(3,3,nthread)
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    logical                      :: shake_end
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: atm1, atm2, iatm(1:8)
    integer                      :: iteration, omp_get_thread_num
    ! FEP
    real(dp)                     :: lambbond

    real(dp),            pointer :: HGr_bond_dist(:,:,:,:)
    real(dp),            pointer :: HGr_bond_vector(:,:,:,:,:)
    real(dp),            pointer :: HGr_shake_force(:,:,:,:)
    real(dp),            pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell


    mass            => domain%mass
    ncell           => domain%num_cell_local

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list 
    HGr_bond_dist   => constraints%HGr_bond_dist
    HGr_bond_vector => constraints%HGr_bond_vector
    HGr_shake_force => constraints%HGr_shake_force

    iteration       =  constraints%shake_iteration
    tolerance       =  real(constraints%shake_tolerance,dp)
    connect         =  constraints%connect

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, atm1, atm2, shake_end, r, x12, y12, z12,  &
    !$omp         dist2, imass1, imass2, x12_old, y12_old, z12_old,            &
    !$omp         factor, g12, g12m1, g12m2, v12m1, v12m2, fx, fy, fz,         &
    !$omp         viri, id, coord_dtmp, vel_dtmp, iatm, lambbond)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! initialize force and store old bond vector
    !
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          atm1 = HGr_bond_list(1,k,j,icel)
          do ih = 1, j
            atm2 = HGr_bond_list(ih+1,k,j,icel)
            HGr_shake_force(ih,k,j,icel) = 0.0_dp
            HGr_bond_vector(1:3,ih,k,j,icel) = &
              real(coord_old(1:3,atm1,icel)-coord_old(1:3,atm2,icel), dp)
          end do
        end do
      end do
    end do

    ! shake iteration
    !
    shake_end = .true.
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / real(mass(iatm(1),icel),dp)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            shake_end = .true.
            do ih = 1, j

              imass2 = 1.0_dp / real(mass(iatm(ih+1),icel),dp)

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12 = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12 = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12 = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              dist2 = x12**2 + y12**2 + z12**2

              if (abs(dist2 - r*r) >= 2.0_dp*tolerance*r) then 

                shake_end = .false.
                x12_old = HGr_bond_vector(1,ih,k,j,icel)
                y12_old = HGr_bond_vector(2,ih,k,j,icel)
                z12_old = HGr_bond_vector(3,ih,k,j,icel)

                factor = (x12*x12_old + y12*y12_old + z12*z12_old)* &
                         (imass1 + imass2)
                g12    = 0.5_dp*(dist2 - r*r)/factor
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2
                v12m1  = g12m1 / real(dt,dp)
                v12m2  = g12m2 / real(dt,dp)

                HGr_shake_force(ih,k,j,icel) = &
                     HGr_shake_force(ih,k,j,icel) + real(g12,dp)

                coord_dtmp(1,1)    = coord_dtmp(1,1)    - g12m1 * x12_old
                coord_dtmp(2,1)    = coord_dtmp(2,1)    - g12m1 * y12_old
                coord_dtmp(3,1)    = coord_dtmp(3,1)    - g12m1 * z12_old
                coord_dtmp(1,ih+1) = coord_dtmp(1,ih+1) + g12m2 * x12_old
                coord_dtmp(2,ih+1) = coord_dtmp(2,ih+1) + g12m2 * y12_old
                coord_dtmp(3,ih+1) = coord_dtmp(3,ih+1) + g12m2 * z12_old

                if (vel_update) then
                  vel_dtmp(1,1)    = vel_dtmp(1,1)    - v12m1 * x12_old
                  vel_dtmp(2,1)    = vel_dtmp(2,1)    - v12m1 * y12_old
                  vel_dtmp(3,1)    = vel_dtmp(3,1)    - v12m1 * z12_old
                  vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + v12m2 * x12_old
                  vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + v12m2 * y12_old
                  vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + v12m2 * z12_old
                end if 

              end if
            end do

            if (shake_end) exit

          end do

          if (.not.shake_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Shake> SHAKE algorithm failed to converge (see "Chapter: Trouble shooting" in the user manual): indexes', &
              domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          endif

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),dp)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),dp)

        end do
      end do
    end do


    ! compute constraint virial
    !   Note: virial =>  virial/dt**2 in compute_constraint
    !
    viri(1:3) = 0.0_dp

    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          do ih = 1, j

            x12_old = HGr_bond_vector(1,ih,k,j,icel)
            y12_old = HGr_bond_vector(2,ih,k,j,icel)
            z12_old = HGr_bond_vector(3,ih,k,j,icel)
            g12 = real(HGr_shake_force(ih,k,j,icel),dp)
      
            fx  = g12 * x12_old
            fy  = g12 * y12_old
            fz  = g12 * z12_old
            viri(1) = viri(1) + x12_old * fx
            viri(2) = viri(2) + y12_old * fy
            viri(3) = viri(3) + z12_old * fz

          end do
        end do
      end do
    end do

    ! FEP: virials for atoms in single topology are scaled by lambbond
    if (domain%fep_use) then
      do icel = id+1, ncell, nthread
        do j = 1, connect
          do k = 1, HGr_local(j,icel)
            iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
            do ih = 1, j
              if ((domain%fepgrp(iatm(ih+1),icel) == 1) .and. &
                  (domain%fepgrp(iatm(1),icel) == 1)) then
                lambbond = real(domain%lambbondA,dp)
              else if ((domain%fepgrp(iatm(ih+1),icel) == 2) .and. &
                       (domain%fepgrp(iatm(1),icel) == 2)) then
                lambbond = real(domain%lambbondB,dp)
              else
                lambbond = 1.0_dp
              end if
              x12_old = HGr_bond_vector(1,ih,k,j,icel)
              y12_old = HGr_bond_vector(2,ih,k,j,icel)
              z12_old = HGr_bond_vector(3,ih,k,j,icel)
              g12 = real(HGr_shake_force(ih,k,j,icel),dp) * (lambbond - 1.0_dp)
              fx  = g12 * x12_old
              fy  = g12 * y12_old
              fz  = g12 * z12_old
              viri(1) = viri(1) + x12_old * fx
              viri(2) = viri(2) + y12_old * fy
              viri(3) = viri(3) + z12_old * fz
            end do
          end do
        end do
      end do
    end if

    virial_omp(1,1,id+1) = virial_omp(1,1,id+1) - viri(1)
    virial_omp(2,2,id+1) = virial_omp(2,2,id+1) - viri(2)
    virial_omp(3,3,id+1) = virial_omp(3,3,id+1) - viri(3)

    !$omp end parallel

    do i = 1, nthread
      virial(1:3,1:3) = virial(1:3,1:3) + real(virial_omp(1:3,1:3,i),dp)
    end do

    return

  end subroutine compute_shake

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_fast_vv1
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_fast_vv1(dt, coord_old, domain, &
                                     constraints, coord, vel)

    ! formal arguments
    real(dp),                 intent(in)    :: dt
    real(dp),                 intent(in)    :: coord_old(:,:,:)
    type(s_domain),   target, intent(in)    :: domain
    type(s_constraints),      intent(inout) :: constraints
    real(dp),                 intent(inout) :: coord(:,:,:)
    real(dp),                 intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                  :: rHH, rOH, massO, massH
    real(dp)                  :: mo, mh, mohh, ra, rb, rc, rc2, tmp, tmp2
    real(dp)                  :: alpa, beta, gama, sqr_ab
    real(dp)                  :: axlng, aylng, azlng
    real(dp)                  :: t1(1:3), t2(1:3), t3(1:3)
    real(dp)                  :: cosphi, costhe, sinphi
    real(dp)                  :: sinthe, cospsi, sinpsi
    real(dp)                  :: aksxd(1:3), aksyd(1:3), akszd(1:3)
    real(dp)                  :: rab(1:3), rac(1:3), rabc(1:3)
    real(dp)                  :: a1(1:3), b1(1:3), c1(1:3)
    real(dp)                  :: a3(1:3), b3(1:3), c3(1:3)
    real(dp)                  :: b0d(1:3), c0d(1:3)
    real(dp)                  :: a1d(1:3), b1d(1:3), c1d(1:3)
    real(dp)                  :: ya2d, yb2d, yc2d
    real(dp)                  :: a3d(1:3), b3d(1:3), c3d(1:3)
    real(dp)                  :: xb2d, xb2d2, hh2, deltx, hhhh
    integer                   :: i, ix, iatom(3)
    integer                   :: id, omp_get_thread_num

    real(dp),         pointer :: mass(:,:)
    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: ncell



    mass       => domain%mass
    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        =  constraints%water_rOH
    rHH        =  constraints%water_rHH
    massO      =  constraints%water_massO
    massH      =  constraints%water_massH

    mo         =  massO
    mh         =  massH
    mohh       =  massO + 2.0_dp * massH
    rc         =  rHH / 2.0_dp
    ra         =  2.0_dp * mh * dsqrt(rOH * rOH - rc * rc)/mohh
    rb         =  dsqrt(rOH * rOH - rc * rc) - ra
    rc2        =  rHH
    hhhh       =  rHH * rHH
    mo         =  mo/mohh
    mh         =  mh/mohh

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, iatom, rab, rac, rabc, a1, b1, c1, a3, b3, c3,        &
    !$omp         aksXd, aksYd, aksZd, axlng, aylng, azlng, t1, t2, t3,        &
    !$omp         b0d, c0d, a1d, b1d, c1d, a3d, b3d, c3d, ya2d, yb2d, yc2d,    &
    !$omp         xb2d, xb2d2, sinphi, cosphi, sinpsi, cospsi, sinthe, costhe, &
    !$omp         tmp, tmp2, hh2, deltx, alpa, beta, gama, sqr_ab, id)
    !
    ! note: reduction cannot be used for "virial" because of the rounding error,
    ! which is very large especialy in the case of langevin NPT.
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)

      ! A1'
      rab(1:3) = coord_old(1:3,iatom(2),i) - coord_old(1:3,iatom(1),i)
      rac(1:3) = coord_old(1:3,iatom(3),i) - coord_old(1:3,iatom(1),i)

      rabc(1:3) = coord(1:3,iatom(1),i)*mo + &
                 (coord(1:3,iatom(2),i) + coord(1:3,iatom(3),i))*mh

      a1(1:3) = coord(1:3,iatom(1),i) - rabc(1:3)
      b1(1:3) = coord(1:3,iatom(2),i) - rabc(1:3)
      c1(1:3) = coord(1:3,iatom(3),i) - rabc(1:3)

      aksZd(1) = rab(2)*rac(3) - rab(3)*rac(2)
      aksZd(2) = rab(3)*rac(1) - rab(1)*rac(3)
      aksZd(3) = rab(1)*rac(2) - rab(2)*rac(1)
      aksXd(1) = a1(2)*aksZd(3) - a1(3)*aksZd(2)
      aksXd(2) = a1(3)*aksZd(1) - a1(1)*aksZd(3)
      aksXd(3) = a1(1)*aksZd(2) - a1(2)*aksZd(1)
      aksYd(1) = aksZd(2)*aksXd(3) - aksZd(3)*aksXd(2)
      aksYd(2) = aksZd(3)*aksXd(1) - aksZd(1)*aksXd(3)
      aksYd(3) = aksZd(1)*aksXd(2) - aksZd(2)*aksXd(1)

      axlng = 1.0_dp/dsqrt(aksXd(1)*aksXd(1)+aksXd(2)*aksXd(2)+aksXd(3)*aksXd(3))
      aylng = 1.0_dp/dsqrt(aksYd(1)*aksYd(1)+aksYd(2)*aksYd(2)+aksYd(3)*aksYd(3))
      azlng = 1.0_dp/dsqrt(aksZd(1)*aksZd(1)+aksZd(2)*aksZd(2)+aksZd(3)*aksZd(3))

      t1(1:3) = aksXd(1:3) * axlng
      t2(1:3) = aksYd(1:3) * aylng
      t3(1:3) = aksZd(1:3) * azlng

      b0d(1) = t1(1)*rab(1) + t1(2)*rab(2) + t1(3)*rab(3)
      b0d(2) = t2(1)*rab(1) + t2(2)*rab(2) + t2(3)*rab(3)
      c0d(1) = t1(1)*rac(1) + t1(2)*rac(2) + t1(3)*rac(3)
      c0d(2) = t2(1)*rac(1) + t2(2)*rac(2) + t2(3)*rac(3)
      a1d(3) = t3(1)*a1(1) + t3(2)*a1(2) + t3(3)*a1(3)
      b1d(1) = t1(1)*b1(1) + t1(2)*b1(2) + t1(3)*b1(3)
      b1d(2) = t2(1)*b1(1) + t2(2)*b1(2) + t2(3)*b1(3)
      b1d(3) = t3(1)*b1(1) + t3(2)*b1(2) + t3(3)*b1(3)
      c1d(1) = t1(1)*c1(1) + t1(2)*c1(2) + t1(3)*c1(3)
      c1d(2) = t2(1)*c1(1) + t2(2)*c1(2) + t2(3)*c1(3)
      c1d(3) = t3(1)*c1(1) + t3(2)*c1(2) + t3(3)*c1(3)

      sinphi = a1d(3) / ra
      tmp    = 1.0_dp - sinphi*sinphi
      if (tmp <= 0.0_dp) then
        cosphi = 0.0_dp
      else
        cosphi = dsqrt(tmp)
      end if

      sinpsi = ( b1d(3) - c1d(3) ) / (rc2 * cosphi)
      tmp2   = 1.0_dp - sinpsi*sinpsi
      if (tmp2 <= 0.0_dp ) then
        cospsi = 0.0_dp
      else
        cospsi = dsqrt(tmp2)
      end if

      ya2d  =   ra * cosphi
      xb2d  = - rc * cospsi
      yb2d  = - rb * cosphi - rc *sinpsi * sinphi
      yc2d  = - rb * cosphi + rc *sinpsi * sinphi
      xb2d2 = xb2d * xb2d
      hh2   = 4.0_dp * xb2d2 + (yb2d-yc2d) * (yb2d-yc2d) &
                           + (b1d(3)-c1d(3)) * (b1d(3)-c1d(3))
      deltx = 2.0_dp * xb2d + dsqrt(4.0_dp * xb2d2 - hh2 + hhhh)
      xb2d  = xb2d - deltx * 0.5_dp


      ! alpha, beta, gamma
      alpa = xb2d * (b0d(1)-c0d(1)) + b0d(2) * yb2d + c0d(2) * yc2d
      beta = xb2d * (c0d(2)-b0d(2)) + b0d(1) * yb2d + c0d(1) * yc2d
      gama = b0d(1)*b1d(2) - b1d(1)*b0d(2) + c0d(1)*c1d(2) - c1d(1)*c0d(2)

      sqr_ab = alpa * alpa + beta * beta
      sinthe = (alpa*gama - beta * dsqrt(sqr_ab - gama * gama)) / sqr_ab


      ! A3'
      costhe = dsqrt(1.0_dp - sinthe * sinthe)
      a3d(1) = - ya2d * sinthe
      a3d(2) =   ya2d * costhe
      a3d(3) = a1d(3)
      b3d(1) =   xb2d * costhe - yb2d * sinthe
      b3d(2) =   xb2d * sinthe + yb2d * costhe
      b3d(3) = b1d(3)
      c3d(1) = - xb2d * costhe - yc2d * sinthe
      c3d(2) = - xb2d * sinthe + yc2d * costhe
      c3d(3) = c1d(3)


      ! A3
      a3(1:3) = t1(1:3)*a3d(1) + t2(1:3)*a3d(2) + t3(1:3)*a3d(3)
      b3(1:3) = t1(1:3)*b3d(1) + t2(1:3)*b3d(2) + t3(1:3)*b3d(3)
      c3(1:3) = t1(1:3)*c3d(1) + t2(1:3)*c3d(2) + t3(1:3)*c3d(3)

      coord(1:3,iatom(1),i) = rabc(1:3) + a3(1:3)
      coord(1:3,iatom(2),i) = rabc(1:3) + b3(1:3)
      coord(1:3,iatom(3),i) = rabc(1:3) + c3(1:3)

      vel(1:3,iatom(1),i) = vel(1:3,iatom(1),i) + (a3(1:3)-a1(1:3))/dt
      vel(1:3,iatom(2),i) = vel(1:3,iatom(2),i) + (b3(1:3)-b1(1:3))/dt
      vel(1:3,iatom(3),i) = vel(1:3,iatom(3),i) + (c3(1:3)-c1(1:3))/dt

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_fast_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_vv1
  !> @brief        RATTLE_VV1 for rigid bonds
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[in]    coord_old   : reference coordinates
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_vv1(dt, coord_old, domain, constraints, coord, vel)

    ! formal arguments
    real(dp),                    intent(in)    :: dt
    real(dp),                    intent(in)    :: coord_old(:,:,:)
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(dp),                    intent(inout) :: coord(:,:,:)
    real(dp),                    intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: x12_old, y12_old, z12_old
    real(dp)                     :: dist, dist2, r
    real(dp)                     :: factor, g12, g12m1, g12m2, v12m1, v12m2
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: atm1, atm2, iatm(1:8)
    integer                      :: iteration, omp_get_thread_num
    logical                      :: shake_end

    real(dp),            pointer :: HGr_bond_vector(:,:,:,:,:)
    real(dp),            pointer :: HGr_bond_dist(:,:,:,:)
    real(dp),            pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell


    mass            => domain%mass
    ncell           => domain%num_cell_local

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    HGr_bond_dist   => constraints%HGr_bond_dist
    HGr_bond_vector => constraints%HGr_bond_vector

    iteration       =  constraints%shake_iteration
    tolerance       =  real(constraints%shake_tolerance,dp)
    connect         =  constraints%connect

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, atm1, atm2, shake_end, r, x12, y12, z12,  &
    !$omp         dist2, dist, imass1, imass2, x12_old, y12_old, z12_old,      &
    !$omp         factor, g12, g12m1, g12m2, v12m1, v12m2,          &
    !$omp         id, coord_dtmp, vel_dtmp, iatm)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! shake iteration
    !
    ! initialize force and store old bond vector
    !
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)
          atm1 = HGr_bond_list(1,k,j,icel)
          do ih = 1, j
            atm2 = HGr_bond_list(ih+1,k,j,icel)
            HGr_bond_vector(1:3,ih,k,j,icel) = &
              real(coord_old(1:3,atm1,icel) - coord_old(1:3,atm2,icel), dp)
          end do
        end do
      end do
    end do

    shake_end = .true.

    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / real(mass(iatm(1),icel),dp)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            shake_end = .true.
            do ih = 1, j

              imass2 = 1.0_dp / real(mass(iatm(ih+1),icel),dp)

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12 = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12 = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12 = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              dist2 = x12**2 + y12**2 + z12**2

              if (abs(dist2 - r*r) >= 2.0_dp*tolerance*r) then

                shake_end = .false.
                x12_old = HGr_bond_vector(1,ih,k,j,icel)
                y12_old = HGr_bond_vector(2,ih,k,j,icel)
                z12_old = HGr_bond_vector(3,ih,k,j,icel)

                factor = (x12*x12_old + y12*y12_old + z12*z12_old)* &
                         (imass1 + imass2)
                g12    = 0.5_dp*(dist2 - r*r)/factor
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2
                v12m1  = g12m1 / real(dt,dp)
                v12m2  = g12m2 / real(dt,dp)

                coord_dtmp(1,1)    = coord_dtmp(1,1)    - g12m1 * x12_old
                coord_dtmp(2,1)    = coord_dtmp(2,1)    - g12m1 * y12_old
                coord_dtmp(3,1)    = coord_dtmp(3,1)    - g12m1 * z12_old
                coord_dtmp(1,ih+1) = coord_dtmp(1,ih+1) + g12m2 * x12_old
                coord_dtmp(2,ih+1) = coord_dtmp(2,ih+1) + g12m2 * y12_old
                coord_dtmp(3,ih+1) = coord_dtmp(3,ih+1) + g12m2 * z12_old

                vel_dtmp(1,1)    = vel_dtmp(1,1)    - v12m1 * x12_old
                vel_dtmp(2,1)    = vel_dtmp(2,1)    - v12m1 * y12_old
                vel_dtmp(3,1)    = vel_dtmp(3,1)    - v12m1 * z12_old
                vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + v12m2 * x12_old
                vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + v12m2 * y12_old
                vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + v12m2 * z12_old

              end if
            end do

            if (shake_end) exit

          end do

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),dp)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),dp)

          if (.not. shake_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Rattle_VV1> SHAKE algorithm failed to converge (see "Chapter: Trouble shooting" in the user manual): '//&
              'indexes',domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          end if

        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_fast_vv2
  !> @brief        SETTLE for three-site water
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_fast_vv2(domain, constraints, coord, vel)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_constraints),      intent(inout) :: constraints
    real(dp),                 intent(inout) :: coord(:,:,:)
    real(dp),                 intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                  :: rHH, rOH
    real(dp)                  :: massO, massH, massOH, massHH, mass2OH
    real(dp)                  :: mo, mh, woh, mohh
    real(dp)                  :: rab(1:3), rbc(1:3), rca(1:3)
    real(dp)                  :: vab(1:3), vbc(1:3), vca(1:3)
    real(dp)                  :: unit_ab(1:3), unit_bc(1:3), unit_ca(1:3)
    real(dp)                  :: length_ab, length_bc, length_ca
    real(dp)                  :: vab0, vbc0, vca0, cosA, cosB, cosC
    real(dp)                  :: MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det
    integer                   :: i, ix, iatom(1:3), IH1, IH2, k
    integer                   :: id, omp_get_thread_num

    integer,          pointer :: nwater(:), water_list(:,:,:)
    integer,          pointer :: ncell


    nwater     => domain%num_water
    water_list => domain%water_list
    ncell      => domain%num_cell_local

    rOH        = constraints%water_rOH
    rHH        = constraints%water_rHH
    massO      = constraints%water_massO
    massH      = constraints%water_massH

    massOH     = massO + massH
    massHH     = massH + massH
    mass2OH    = 2.0_dp * massOH

    !$omp parallel default(shared)                                             &
    !$omp private(i, ix, k, rab, rbc, rca, vab, vbc, vca, length_ab, length_bc,&
    !$omp         length_ca, unit_ab, unit_bc, unit_ca, vab0, vbc0, vca0,      &
    !$omp         cosA, cosB, cosC, MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det,   &
    !$omp         iatom, id)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, nwater(i)

      iatom(1:3) = water_list(1:3,ix,i)

      rab(1:3) = coord(1:3,iatom(2),i) - coord(1:3,iatom(1),i)
      rbc(1:3) = coord(1:3,iatom(3),i) - coord(1:3,iatom(2),i)
      rca(1:3) = coord(1:3,iatom(1),i)  - coord(1:3,iatom(3),i)

      vab(1:3) = vel(1:3,iatom(2),i) - vel(1:3,iatom(1),i)
      vbc(1:3) = vel(1:3,iatom(3),i) - vel(1:3,iatom(2),i)
      vca(1:3) = vel(1:3,iatom(1),i)  - vel(1:3,iatom(3),i)

      length_ab = dsqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)
      length_bc = dsqrt(rbc(1)**2 + rbc(2)**2 + rbc(3)**2)
      length_ca = dsqrt(rca(1)**2 + rca(2)**2 + rca(3)**2)

      unit_ab(1:3) = rab(1:3)/length_ab
      unit_bc(1:3) = rbc(1:3)/length_bc
      unit_ca(1:3) = rca(1:3)/length_ca

      vab0 = 0.0_dp
      vbc0 = 0.0_dp
      vca0 = 0.0_dp

      do k = 1, 3
        vab0 = vab0 + vab(k)*unit_ab(k)
        vbc0 = vbc0 + vbc(k)*unit_bc(k)
        vca0 = vca0 + vca(k)*unit_ca(k)
      end do

      cosA = 0.0_dp
      cosB = 0.0_dp
      cosC = 0.0_dp

      do k = 1, 3
        cosA = cosA - unit_ab(k)*unit_ca(k)
        cosB = cosB - unit_bc(k)*unit_ab(k)
        cosC = cosC - unit_ca(k)*unit_bc(k)
      end do

      MbCA = massH*cosC*cosA - massOH*cosB
      MaBC = massO*cosB*cosC - massHH*cosA
      McAB = massH*cosA*cosB - massOH*cosC

      T_AB = vab0*(mass2OH-massO*cosC**2) + vbc0*MbCA + vca0*MaBC
      T_AB = T_AB * massO
      T_BC = vbc0*(massOH**2-(massH*cosA)**2)
      T_BC = T_BC + vca0*massO*McAB + vab0*massO*MbCA
      T_CA = vca0*(mass2OH-massO*cosB**2) + vab0*MaBC + vbc0*McAB
      T_CA = T_CA * massO
      Det  = 2.0_dp*massOH**2 + 2.0_dp*massO*massH*cosA*cosB*cosC
      Det  = Det - 2.0_dp*massH**2*cosA**2 - massO*massOH*(cosB**2+cosC**2)
      Det  = Det / massH

      T_AB = T_AB / Det
      T_BC = T_BC / Det
      T_CA = T_CA / Det

      vel(1:3,iatom(1),i)  = vel(1:3,iatom(1),i) +  &
                       (T_AB*unit_ab(1:3)-T_CA*unit_ca(1:3))/massO
      vel(1:3,iatom(2),i) = vel(1:3,iatom(2),i) + &
                       (T_BC*unit_bc(1:3)-T_AB*unit_ab(1:3))/massH
      vel(1:3,iatom(3),i) = vel(1:3,iatom(3),i) + &
                       (T_CA*unit_ca(1:3)-T_BC*unit_bc(1:3))/massH

      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_fast_vv2
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_rattle_vv2
  !> @brief        RATTLE_VV2 for rigid bonds
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] coord       : coordinates
  !! @param[inout] vel         : velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_rattle_vv2(domain, constraints, coord, vel)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    real(dp),                    intent(inout) :: coord(:,:,:)
    real(dp),                    intent(inout) :: vel(:,:,:)

    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: vel_x12, vel_y12, vel_z12
    real(dp)                     :: r, dot_d12_v12
    real(dp)                     :: factor, g12, g12m1, g12m2
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: iteration, omp_get_thread_num
    integer                      :: iatm(1:8)
    logical                      :: rattle_end

    real(dp),            pointer :: HGr_bond_dist(:,:,:,:)
    real(dp),            pointer :: mass(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer,             pointer :: ncell


    mass          => domain%mass
    ncell         => domain%num_cell_local

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    iteration     =  constraints%shake_iteration
    tolerance     =  real(constraints%shake_tolerance,dp)
    connect       =  constraints%connect

    !$omp parallel default(shared)                                             &
    !$omp private(icel, i, j, k, ih, rattle_end, r, x12, y12, z12,             &
    !$omp         dot_d12_v12, imass1, imass2, vel_x12, vel_y12, vel_z12,      &
    !$omp         factor, g12, g12m1, g12m2, id, coord_dtmp, vel_dtmp, iatm)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! rattle iteration
    !
    do icel = id+1, ncell, nthread
      do j = 1, connect
        do k = 1, HGr_local(j,icel)

          iatm(1:j+1) = HGr_bond_list(1:j+1,k,j,icel)
          imass1 = 1.0_dp / mass(iatm(1),icel)

          coord_dtmp(1:3,1:j+1) = &
            real(coord(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(vel  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            rattle_end = .true.

            do ih = 1, j

              imass2 = 1.0_dp / dble(mass(iatm(ih+1),icel))

              r = real(HGr_bond_dist(ih+1,k,j,icel),dp)

              x12     = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12     = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12     = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              vel_x12 = vel_dtmp(1,1)   - vel_dtmp(1,ih+1)
              vel_y12 = vel_dtmp(2,1)   - vel_dtmp(2,ih+1)
              vel_z12 = vel_dtmp(3,1)   - vel_dtmp(3,ih+1)

              dot_d12_v12 = x12*vel_x12 + y12*vel_y12 + z12*vel_z12

              if (abs(dot_d12_v12) >= tolerance) then

                rattle_end = .false.

                g12 = dot_d12_v12 / ( (imass1+imass2)*r**2)
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2

                vel_dtmp(1,1)    = vel_dtmp(1,1)    - g12m1 * x12
                vel_dtmp(2,1)    = vel_dtmp(2,1)    - g12m1 * y12
                vel_dtmp(3,1)    = vel_dtmp(3,1)    - g12m1 * z12
                vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + g12m2 * x12
                vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + g12m2 * y12
                vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + g12m2 * z12

              end if
            end do

            if (rattle_end) exit

          end do

          coord(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),dp)

          vel  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),dp)

          if (.not. rattle_end) then
            write(MsgOut, '(A,8i10)') &
              'Compute_Rattle_VV2> SHAKE algorithm failed to converge (see "Chapter: Trouble shooting" in the user manual): '//&
              'indexes',domain%id_l2g(iatm(1:j+1), icel)
            call error_msg('')
          end if

        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine compute_rattle_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    decide_dummy
  !> @brief        Decide Dummy atom posiion
  !! @authors      JJ
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine decide_dummy(domain, constraints, coord)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints),         intent(in)    :: constraints
    real(dp),                    intent(inout) :: coord(:,:,:)

    integer,       pointer :: ncell, nwater(:), water_list(:,:,:)

    ! local variables
    integer                :: i, j, index(4)
    integer                :: id, omp_get_thread_num
    real(dp)               :: rOD, dist_rijk, rij(3), rjk(3)
    real(dp)               :: rijk(3), uijk(3)

    ncell          => domain%num_cell_local
    nwater         => domain%num_water
    water_list     => domain%water_list

    rOD        =  constraints%water_rOD

    !$omp parallel default(shared)                                     &
    !$omp          private(id, i, j, index, rij, rjk, rijk, dist_rijk, &
    !$omp                  uijk)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do j = 1, nwater(i)

        index(1:4) = water_list(1:4,j,i)
        rij(1:3)   = coord(1:3,index(2),i) - coord(1:3,index(1),i)
        rjk(1:3)   = coord(1:3,index(3),i) - coord(1:3,index(2),i)
        rijk(1:3)  = rij(1:3) + 0.5_wp*rjk(1:3)
        dist_rijk  = sqrt(rijk(1)*rijk(1) + rijk(2)*rijk(2) + rijk(3)*rijk(3))
        uijk(1:3)  = rijk(1:3) / dist_rijk
        coord(1:3,index(4),i) = coord(1:3,index(1),i) &
                              + rOD*uijk(1:3)
      end do
    end do
    !$omp end parallel

    return

  end subroutine decide_dummy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    water_force_redistribution
  !> @brief        Dummy atom force is redistributed to other OH2 atoms
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine water_force_redistribution(constraints, domain, force, virial)

    ! formal arguments
    type(s_constraints),         intent(in)    :: constraints
    type(s_domain),      target, intent(in)    :: domain
    real(dp),                    intent(inout) :: force(:,:,:)
    real(dp),                    intent(inout) :: virial(:,:)

    integer,       pointer :: ncell, nwater(:), water_list(:,:,:)
    real(dp),      pointer :: coord(:,:,:)

    ! local variables
    integer                :: i, j, index(4), k
    integer                :: id, omp_get_thread_num
    real(dp)               :: factor, rOD, dist_rijk, dist_rid
    real(dp)               :: rid(3), rij(3), rjk(3), rijk(3)
    real(dp)               :: Fd(3), F1(3), rid_Fd
    real(dp)               :: virial_sub(3), virial_add(3)

    ncell          => domain%num_cell_local
    nwater         => domain%num_water
    water_list     => domain%water_list
    coord          => domain%coord

    rOD        =  constraints%water_rOD

    virial_sub(1:3) = 0.0_dp
    virial_add(1:3) = 0.0_dp

    !$omp parallel default(shared)                                          &
    !$omp          private(id, i, j, index, rid, rij, rjk, rijk, dist_rijk, &
    !$omp                  dist_rid, factor, Fd, rid_Fd, F1)                &
    !$omp          reduction(+:virial_sub, virial_add)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread

      do j = 1, nwater(i)
        index(1:4) = water_list(1:4,j,i)
        do k = 1, 4
          virial_sub(1:3) = virial_sub(1:3) &
                          + force(1:3,index(k),i)*coord(1:3,index(k),i)
        end do
        rid(1:3)   = coord(1:3,index(4),i) - coord(1:3,index(1),i)
        rij(1:3)   = coord(1:3,index(2),i) - coord(1:3,index(1),i)
        rjk(1:3)   = coord(1:3,index(3),i) - coord(1:3,index(2),i)
        rijk(1:3)  = rij(1:3) + 0.5_wp*rjk(1:3)
        dist_rijk  = sqrt(rijk(1)*rijk(1) + rijk(2)*rijk(2) + rijk(3)*rijk(3))
        dist_rid   = rid(1)*rid(1) + rid(2)*rid(2) + rid(3)*rid(3)

        factor     = rOD / dist_rijk
        Fd(1:3)    = force(1:3,index(4),i)
        rid_Fd     = rid(1)*Fd(1) + rid(2)*Fd(2) + rid(3)*Fd(3)
        F1(1:3)    = rid(1:3) * (rid_Fd/dist_rid)
        F1(1:3)    = Fd(1:3) - F1(1:3)
        force(1:3,index(1),i) = force(1:3,index(1),i) &
                              + Fd(1:3) - factor*F1(1:3)
        force(1:3,index(2),i) = force(1:3,index(2),i) &
                              + 0.5_wp*factor*F1(1:3)
        force(1:3,index(3),i) = force(1:3,index(3),i) &
                              + 0.5_wp*factor*F1(1:3)
        do k = 1, 3
          virial_add(1:3) = virial_add(1:3) &
                          + force(1:3,index(k),i)*coord(1:3,index(k),i)
        end do
     end do

   end do
   !$omp end parallel

   virial(1,1) = virial(1,1) - virial_sub(1) + virial_add(1)
   virial(2,2) = virial(2,2) - virial_sub(2) + virial_add(2)
   virial(3,3) = virial(3,3) - virial_sub(3) + virial_add(3)

   end subroutine water_force_redistribution

end module sp_constraints_mod
